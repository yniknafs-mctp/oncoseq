'''
Created on Aug 7, 2011

@author: mkiyer
'''
import argparse
import logging
import sys
import os

from oncoseq.lib import config
from oncoseq.lib.config import AnalysisConfig, PipelineConfig
from oncoseq.lib.cluster import submit_job_pbs, submit_job_nopbs
from oncoseq.lib.base import up_to_date
from oncoseq.lib.defuse import get_defuse_config_string

# setup pipeline script files
import oncoseq.rnaseq.pipeline
_pipeline_dir = oncoseq.rnaseq.pipeline.__path__[0] 

def run_lane(lane, genome, server, pipeline, num_processors,
             submit_job_func):
    #
    # create tmp and log directories
    #
    tmp_dir = os.path.join(lane.output_dir, "tmp") 
    if not os.path.exists(tmp_dir):
        logging.info("Creating directory: %s" % (tmp_dir))
        os.makedirs(tmp_dir)
    log_dir = os.path.join(lane.output_dir, "log")
    if not os.path.exists(log_dir):
        logging.info("Creating directory: %s" % (log_dir))
        os.makedirs(log_dir)
    #
    # copy and uncompress reads
    # 
    msg = "Copying and uncompressing read sequences"
    copy_fastq_deps = []
    for readnum in xrange(len(lane.fastq_files)):
        if up_to_date(lane.copied_fastq_files[readnum], lane.fastq_files[readnum]):
            logging.info("[SKIPPED] %s read%d" % (msg, readnum+1))
        else:
            logging.info("%s read%d" % (msg, readnum+1))       
            pyscript = os.path.join(_pipeline_dir, "copy_uncompress_fastq.py")    
            args = [sys.executable, pyscript,
                    "--quals", lane.quality_scores,
                    lane.fastq_files[readnum],
                    lane.copied_fastq_files[readnum]]
            log_file = os.path.join(log_dir, "copy_uncompress_fastq_read%d.log" % (readnum+1))
            job_id = submit_job_func("cp%d_%s" % (readnum+1, lane.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=lane.output_dir,
                                     walltime="20:00:00",
                                     stderr_filename=log_file)
            copy_fastq_deps.append(job_id)
    #
    # run fastqc program
    #
    skip_run = True
    for readnum in xrange(len(lane.fastq_files)):    
        skip_run = skip_run and (up_to_date(lane.fastqc_data_files[readnum], lane.copied_fastq_files[readnum]) and
                                 up_to_date(lane.fastqc_report_files[readnum], lane.copied_fastq_files[readnum]))
    msg = "Running FASTQC quality assessment"
    if skip_run:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info("%s" % msg)       
        num_threads = min(num_processors, len(lane.copied_fastq_files))
        args = [pipeline.fastqc_bin, 
                "--threads", num_threads,
                "-o", lane.output_dir]
        args.extend(lane.copied_fastq_files)        
        log_file = os.path.join(log_dir, "fastqc.log")
        job_id = submit_job_func("fqc_%s" % (lane.id), args,
                                 num_processors=num_threads,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="10:00:00",
                                 deps=copy_fastq_deps,
                                 stderr_filename=log_file)
    #
    # Map reads against abundant sequences
    #    
    msg = "Mapping reads against abundant sequences using bowtie2"
    abundant_mapping_deps = []
    for readnum in xrange(len(lane.copied_fastq_files)):
        abundant_sam_file = lane.abundant_sam_files[readnum]
        if up_to_date(abundant_sam_file, lane.copied_fastq_files[readnum]):
            logging.info("[SKIPPED] %s read%d" % (msg, readnum+1))
        else:
            logging.info("%s read%d" % (msg, readnum+1))
            args = [pipeline.bowtie2_bin, "-p", num_processors, "--phred33",
                    "--end-to-end", "--sensitive", "--reorder", "-M", 200,
                    "-x", os.path.join(server.references_dir, genome.get_path("abundant_bowtie_index")),
                    "-U", lane.copied_fastq_files[readnum],
                    "-S", abundant_sam_file]
            log_file = os.path.join(log_dir, "bowtie2_abundant_seq_read%d.log" % (readnum+1))
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("ab%d_%s" % (readnum+1, lane.id), args,
                                     num_processors=num_processors,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=lane.output_dir,
                                     walltime="20:00:00",
                                     deps=copy_fastq_deps,
                                     stderr_filename=log_file)
            abundant_mapping_deps.append(job_id)
    #
    # filter reads mapping to abundant sequences and output sequences 
    # for tophat alignment
    #
    msg = "Filtering abundant sequences"
    filtered_fastq_deps = []
    skip_run = (up_to_date(lane.abundant_bam_file, f) for f in lane.abundant_sam_files)    
    for readnum in xrange(len(lane.fastq_files)):
        skip_run = skip_run and up_to_date(lane.filtered_fastq_files[readnum], lane.abundant_sam_files[readnum])
    if skip_run:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info("%s" % (msg))
        args = [sys.executable, os.path.join(_pipeline_dir, "filter_abundant_sequences.py"),
                lane.abundant_bam_file]
        for readnum in xrange(len(lane.fastq_files)):
            args.append(lane.abundant_sam_files[readnum])
            args.append(lane.filtered_fastq_files[readnum])
        log_file = os.path.join(log_dir, "filter_abundant_sequences.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("ab_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="10:00:00",
                                 deps=abundant_mapping_deps,
                                 stderr_filename=log_file)
        filtered_fastq_deps = [job_id]

    #
    # run defuse gene fusion prediction
    #
    msg = "Running DeFuse gene fusion prediction"
    if all(up_to_date(lane.defuse_results_file,f) for f in lane.filtered_fastq_files):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info("%s" % (msg))
        if not os.path.exists(lane.defuse_dir):
            logging.info("\tcreating directory: %s" % (lane.defuse_dir))
            os.makedirs(lane.defuse_dir)
        defuse_config = pipeline.defuse_config
        species_root_dir = os.path.join(server.references_dir, genome.root_dir) 
        # write defuse config file
        config_str = get_defuse_config_string(species_root_dir, 
                                              defuse_config,
                                              max_fragment_size=pipeline.max_fragment_size,
                                              num_processors=num_processors)
        f = open(lane.defuse_config_file, "w")
        f.write(config_str)
        f.close()
        # setup command line
        script = os.path.join(defuse_config.source_dir, "scripts", "defuse.pl")
        args = [script, "-c", lane.defuse_config_file, "-d", lane.output_dir,
                "-o", lane.defuse_dir, "-p", num_processors] 
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        log_stderr_file = os.path.join(log_dir, "defuse_stderr.log")
        log_stdout_file = os.path.join(log_dir, "defuse_stdout.log")
        # allocate 12gb to run defuse
        defuse_pmem = int(round(float(12000.0 / num_processors),0))
        job_id = submit_job_func("defuse_%s" % (lane.id), args,
                                 num_processors=num_processors,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 pmem=defuse_pmem,
                                 walltime="80:00:00",
                                 deps=filtered_fastq_deps,
                                 stdout_filename=log_stdout_file,
                                 stderr_filename=log_stderr_file)
    #
    # sort abundant reads bam file
    #
    msg = "Sorting abundant hits BAM file"
    if up_to_date(lane.sorted_abundant_bam_file, lane.abundant_bam_file):    
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "SortSam.jar"),
                "INPUT=%s" % (lane.abundant_bam_file),
                "OUTPUT=%s" % (lane.sorted_abundant_bam_file),
                "SO=coordinate",
                "CREATE_INDEX=true",
                "TMP_DIR=%s" % tmp_dir]
        log_file = os.path.join(log_dir, "picard_sort_abundant.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("sortsam_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 pmem=8192,
                                 walltime="10:00:00",
                                 deps=filtered_fastq_deps,
                                 stderr_filename=log_file)
    #
    # Map reads against contaminant sequences
    #    
    msg = "Mapping reads against foreign contaminant sequences using bowtie2"
    xeno_mapping_deps = []
    for readnum in xrange(len(lane.filtered_fastq_files)):
        xeno_sam_file = lane.xeno_sam_files[readnum]
        if up_to_date(xeno_sam_file, lane.filtered_fastq_files[readnum]):
            logging.info("[SKIPPED] %s read%d" % (msg, readnum+1))
        else:
            logging.info("%s read%d" % (msg, readnum+1))
            args = [pipeline.bowtie2_bin, "-p", num_processors, "--phred33",
                    "--end-to-end", "--sensitive", "--reorder", "-M", 200,
                    "-x", os.path.join(server.references_dir, genome.get_path("xeno_bowtie_index")),
                    "-U", lane.filtered_fastq_files[readnum],
                    "-S", xeno_sam_file]
            log_file = os.path.join(log_dir, "bowtie2_xeno_seq_read%d.log" % (readnum+1))
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("xn%d_%s" % (readnum+1, lane.id), args,
                                     num_processors=num_processors,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=lane.output_dir,
                                     walltime="20:00:00",
                                     deps=filtered_fastq_deps,
                                     stderr_filename=log_file)
            xeno_mapping_deps.append(job_id)
    #
    # filter reads mapping to abundant sequences and output sequences 
    # for tophat alignment
    #
    msg = "Filtering foreign contaminant sequences"
    filter_xeno_deps = []
    if all(up_to_date(lane.xeno_bam_file, f) for f in lane.xeno_sam_files):    
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info("%s" % (msg))
        args = [sys.executable, os.path.join(_pipeline_dir, "filter_unmapped_pairs.py"),
                lane.xeno_bam_file]
        args.extend(lane.xeno_sam_files)
        log_file = os.path.join(log_dir, "filter_xeno_sequences.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("fxn_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="10:00:00",
                                 deps=xeno_mapping_deps,
                                 stderr_filename=log_file)
        filter_xeno_deps = [job_id]
    #
    # sort foreign contaminants bam file
    #
    msg = "Sorting foreign contaminants BAM file"
    if up_to_date(lane.sorted_xeno_bam_file, lane.xeno_bam_file):    
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "SortSam.jar"),
                "INPUT=%s" % (lane.xeno_bam_file),
                "OUTPUT=%s" % (lane.sorted_xeno_bam_file),
                "SO=coordinate",
                "CREATE_INDEX=true",
                "TMP_DIR=%s" % tmp_dir]
        log_file = os.path.join(log_dir, "picard_sort_xeno.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("sortxeno_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="10:00:00",
                                 deps=filter_xeno_deps,
                                 stderr_filename=log_file)
    #
    # determine the fragment size distribution of the reads
    #
    msg = "Estimating fragment size distribution"
    frag_size_deps = []
    if (all(up_to_date(lane.frag_size_dist_file, f) for f in lane.filtered_fastq_files) and
        all(up_to_date(lane.frag_size_dist_plot_file, f) for f in lane.filtered_fastq_files)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info("%s" % (msg))
        args = [sys.executable, 
                os.path.join(_pipeline_dir, "estimate_fragment_size_distribution.py"),
                "--min-fragment-size", pipeline.min_fragment_size,
                "--max-fragment-size", pipeline.max_fragment_size,
                "--bowtie-bin", pipeline.bowtie_bin,
                # TODO: set this?
                "--trim5", 0,
                "--library-fragment-size", lane.library.fragment_length,
                "--default-mean", pipeline.fragment_size_mean_default,
                "--default-stdev", pipeline.fragment_size_stdev_default,
                "--default-adaptor-length", pipeline.adaptor_length_default,
                "-p", num_processors,
                "-n", 1000000,
                os.path.join(server.references_dir, genome.get_path("fragment_size_bowtie_index")),
                lane.frag_size_dist_file,
                lane.frag_size_dist_plot_file]
        args.extend(lane.filtered_fastq_files)   
        log_file = os.path.join(log_dir, "estimate_fragment_size_distribution.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("fs_%s" % (lane.id), args,
                                 num_processors=num_processors,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="10:00:00",
                                 deps=filtered_fastq_deps,
                                 stderr_filename=log_file)
        frag_size_deps = [job_id]
    #
    # align reads with tophat
    #
    tophat_deps = []
    msg = "Aligning reads with Tophat"
    if (up_to_date(lane.tophat_bam_file, lane.frag_size_dist_file) and
        all(up_to_date(lane.tophat_bam_file, f) for f in lane.filtered_fastq_files)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info("%s" % (msg))
        args = [sys.executable, os.path.join(_pipeline_dir, "run_tophat.py"),
                "-p", num_processors,
                "--library-type", config.get_tophat_library_type(lane.library.strand_protocol),
                '--rg-id', lane.id,
                '--rg-sample="%s"' % lane.library.sample.id,
                '--rg-library="%s"' % lane.library.id,
                '--rg-description="%s"' % lane.library.description,
                '--rg-platform-unit="%s"' % lane.lane,
                '--rg-center="%s"' % lane.center_name,
                '--rg-platform="%s"' % lane.platform]
        for arg in pipeline.tophat_args:
            # substitute species-specific root directory
            species_arg = arg.replace("${SPECIES}", os.path.join(server.references_dir, genome.root_dir)) 
            args.extend(['--tophat-arg="%s"' % species_arg])
        args.extend([lane.tophat_dir,
                     os.path.join(server.references_dir, genome.get_path("genome_bowtie_index")),
                     lane.frag_size_dist_file])
        args.extend(lane.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        log_file = os.path.join(log_dir, "tophat.log")
        # allocate 10gb to run tophat
        tophat_pmem = int(round(float(10000.0 / num_processors),0))
        job_id = submit_job_func("tophat_%s" % (lane.id), args,
                                 num_processors=num_processors,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 pmem=tophat_pmem,
                                 walltime="80:00:00",
                                 deps=frag_size_deps,
                                 stderr_filename=log_file)
        tophat_deps = [job_id]
    #
    # run picard diagnostics for alignment results
    #
    msg = "Collecting alignment metrics with Picard"    
    if (up_to_date(lane.alignment_summary_metrics, lane.tophat_bam_file) and
        up_to_date(lane.quality_by_cycle_metrics, lane.tophat_bam_file) and
        up_to_date(lane.quality_distribution_metrics, lane.tophat_bam_file)):            
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "CollectMultipleMetrics.jar"),
                "INPUT=%s" % (lane.tophat_bam_file),
                "REFERENCE_SEQUENCE=%s" % os.path.join(server.references_dir, genome.get_path("genome_fasta_file")),
                "OUTPUT=lane",
                "ASSUME_SORTED=TRUE",
                "TMP_DIR=%s" % tmp_dir,
                "VALIDATION_STRINGENCY=SILENT"]
        log_file = os.path.join(log_dir, "picard_metrics.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        submit_job_func("picard_%s" % (lane.id), args,
                        num_processors=1,
                        node_processors=server.node_processors,
                        node_memory=server.node_mem,
                        pbs_script_lines=server.pbs_script_lines,
                        working_dir=lane.output_dir,
                        pmem=8192,
                        walltime="20:00:00",
                        deps=tophat_deps,
                        stderr_filename=log_file)
    #
    # run picard for rna-seq diagnostics
    #
    msg = "Collecting RNA-Seq metrics with Picard"    
    if (up_to_date(lane.rnaseq_metrics, lane.tophat_bam_file)):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "CollectRnaSeqMetrics.jar"),
                "INPUT=%s" % (lane.tophat_bam_file),
                "REF_FLAT=%s" % (os.path.join(server.references_dir, genome.get_path("gene_annotation_refflat"))),
                "RIBOSOMAL_INTERVALS=%s" % (os.path.join(server.references_dir, genome.get_path("picard_ribosomal_intervals"))),
                "STRAND_SPECIFICITY=%s" % config.get_picard_strand_specificity(lane.library.strand_protocol),                
                "REFERENCE_SEQUENCE=%s" % (os.path.join(server.references_dir, genome.get_path("genome_fasta_file"))),
                "OUTPUT=%s" % (lane.rnaseq_metrics),
                "CHART_OUTPUT=%s" % (lane.rnaseq_metrics_pdf),
                "TMP_DIR=%s" % tmp_dir]
                #"ASSUME_SORTED=TRUE",
                #"VALIDATION_STRINGENCY=SILENT"]
        log_file = os.path.join(log_dir, "picard_rnaseq_metrics.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        submit_job_func("picard_%s" % (lane.id), args,
                        num_processors=1,
                        node_processors=server.node_processors,
                        node_memory=server.node_mem,
                        pbs_script_lines=server.pbs_script_lines,
                        working_dir=lane.output_dir,
                        pmem=8192,
                        walltime="20:00:00",
                        deps=tophat_deps,
                        stderr_filename=log_file)
    #
    # generate genome coverage bedgraph  
    #
    msg = "Generating coverage bedgraph file"
    bedgraph_deps = []
    if (up_to_date(lane.coverage_bedgraph_file, lane.tophat_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [os.path.join(pipeline.bedtools_dir, "genomeCoverageBed"),
                "-bg", "-split", "-g", (os.path.join(server.references_dir, genome.get_path("chrom_sizes"))),
                "-ibam", lane.tophat_bam_file]
        log_file = os.path.join(log_dir, "genomeCoverageBed.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("bedtools_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="20:00:00",
                                 deps=tophat_deps,
                                 stdout_filename=lane.coverage_bedgraph_file,
                                 stderr_filename=log_file)
        bedgraph_deps = [job_id]
    #
    # convert bedgraph to bigwig coverage file
    #
    msg = "Create bigWig file to display coverage"
    if (up_to_date(lane.coverage_bigwig_file, lane.coverage_bedgraph_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [os.path.join(pipeline.ucsc_dir, "bedGraphToBigWig"),
                lane.coverage_bedgraph_file,
                os.path.join(server.references_dir, genome.get_path("chrom_sizes")),
                lane.coverage_bigwig_file]
        log_file = os.path.join(log_dir, "bedGraphToBigWig.log")
        job_id = submit_job_func("bigwig_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="20:00:00",
                                 deps=bedgraph_deps,
                                 stderr_filename=log_file)
    return tophat_deps


def run_library(library, genome, server, pipeline, num_processors,
                submit_job_func):
    lane_deps = []
    has_paired_end = False
    frag_size_dist_files = []
    for lane in library.lanes:
        logging.info("Analyzing lane: %s" % (lane.id))
        #
        # create lane directory
        # 
        if not os.path.exists(lane.output_dir):
            logging.info("Creating directory: %s" % (lane.output_dir))
            os.makedirs(lane.output_dir)
        #
        # process lane
        #
        lane_deps.extend(run_lane(lane, genome, server, pipeline, num_processors, submit_job_func))
        #
        # keep track of whether we have paired end
        #
        frag_size_dist_files.append(lane.frag_size_dist_file)
        if lane.fragment_layout == "paired":
            has_paired_end = True
    # 
    # create directories for tmp and log files
    #
    tmp_dir = os.path.join(library.output_dir, "tmp") 
    if not os.path.exists(tmp_dir):
        logging.info("Creating directory: %s" % (tmp_dir))
        os.makedirs(tmp_dir)
    log_dir = os.path.join(library.output_dir, "log")
    if not os.path.exists(log_dir):
        logging.info("Creating directory: %s" % (log_dir))
        os.makedirs(log_dir)
    #
    # merge lane fragment size distributions
    #
    msg = "Merging lane fragment size distributions"
    merge_frag_size_deps = []
    if all(up_to_date(library.merged_frag_size_dist_file, lane.tophat_bam_file) 
           for lane in library.lanes):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, 
                os.path.join(_pipeline_dir, "merge_fragment_size_distributions.py"),
                library.merged_frag_size_dist_file]
        args.extend(frag_size_dist_files)
        log_file = os.path.join(log_dir, "merge_fragment_size_distributions.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("mfs_%s" % (library.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="10:00:00",
                                 deps=lane_deps,
                                 stderr_filename=log_file)
        merge_frag_size_deps = [job_id]
    #
    # merge lane alignment results
    #
    merge_bam_deps = []
    msg = "Merging lane BAM files"
    if all(up_to_date(library.merged_bam_file, lane.tophat_bam_file) for lane in library.lanes):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "MergeSamFiles.jar")]
        for lane in library.lanes:
            args.append("INPUT=%s" % (lane.tophat_bam_file))
        args.extend(["OUTPUT=%s" % (library.merged_bam_file),
                     "SO=coordinate",
                     "MERGE_SEQUENCE_DICTIONARIES=true",
                     "USE_THREADING=true",
                     "CREATE_INDEX=true",
                     "TMP_DIR=%s" % tmp_dir])
        log_file = os.path.join(log_dir, "picard_merge_bam_files.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("merge_%s" % (library.id), args,
                                 num_processors=min(2, num_processors),
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="20:00:00",
                                 deps=lane_deps,
                                 stderr_filename=log_file)
        merge_bam_deps = [job_id]
    #
    # call snps
    #
    msg = "Calling SNVs with samtools"
    samtools_snv_deps = []
    if up_to_date(library.samtools_vcf_file, library.merged_bam_file):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_pipeline_dir, "call_snps.py"),
                os.path.join(server.references_dir, genome.get_path("genome_fasta_file")),
                library.merged_bam_file,
                library.samtools_bcf_file,
                library.samtools_vcf_file]
        log_file = os.path.join(log_dir, "samtools_snp_calling.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("samsnv_%s" % (library.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="60:00:00",
                                 deps=merge_bam_deps,
                                 stderr_filename=log_file)
        samtools_snv_deps = [job_id]    
    # 
    # Call SNVs using varscan
    #
    msg = "Calling SNVs with VarScan"
    varscan_deps = []
    if up_to_date(library.varscan_snv_file, library.merged_bam_file):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_pipeline_dir, "call_snps_varscan.py"),
                "--varscan-dir", pipeline.varscan_dir,
                os.path.join(server.references_dir, genome.get_path("genome_fasta_file")),
                library.merged_bam_file,
                library.varscan_snv_file,
                library.varscan_indel_file]
        log_stdout_file = os.path.join(log_dir, "varscan_snp_calling_stdout.log")
        log_stderr_file = os.path.join(log_dir, "varscan_snp_calling_stderr.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("varscan_%s" % (library.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="60:00:00",
                                 deps=merge_bam_deps,
                                 stdout_filename=log_stdout_file,
                                 stderr_filename=log_stderr_file)
        varscan_deps = [job_id]
    #
    # run cufflinks to estimate transcript abundance of known genes
    #
    library.cufflinks_gtf_file = os.path.join(library.cufflinks_dir, "transcripts.gtf")    
    cufflinks_deps = []
    msg = "Estimating known transcript abundances with Cufflinks"
    if (up_to_date(library.cufflinks_gtf_file, library.merged_bam_file) and
        up_to_date(library.cufflinks_gtf_file, library.merged_frag_size_dist_file)):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        if not os.path.exists(library.cufflinks_dir):
            logging.info("\tcreating directory: %s" % (library.cufflinks_dir))
            os.makedirs(library.cufflinks_dir)
        args = [sys.executable, os.path.join(_pipeline_dir, "run_cufflinks.py"),
                "--cufflinks-bin", pipeline.cufflinks_bin,
                "-p", num_processors,
                "-L", library.id,
                "--library-type", config.get_tophat_library_type(library.strand_protocol)]
        if has_paired_end:
            args.append("--learn-frag-size")
        for arg in pipeline.cufflinks_args:
            # substitute species-specific root directory
            species_arg = arg.replace("${SPECIES}", os.path.join(server.references_dir, genome.root_dir)),
            args.append('--cufflinks-arg="%s"' % (species_arg))
        args.extend([library.merged_bam_file,
                     library.cufflinks_dir,
                     library.merged_frag_size_dist_file])
        log_file = os.path.join(log_dir, "cufflinks.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        cuff_pmem = int(round(float(24000.0 / num_processors),0))
        job_id = submit_job_func("cuff_%s" % (library.id), args,
                                 num_processors=num_processors,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.cufflinks_dir,
                                 walltime="60:00:00",
                                 pmem=cuff_pmem,
                                 deps=merge_bam_deps + merge_frag_size_deps,
                                 stderr_filename=log_file)
        cufflinks_deps = [job_id]
    return cufflinks_deps + samtools_snv_deps + varscan_deps

def run_sample(sample, genome, server, pipeline, num_processors,
               submit_job_func):
    deps = []
    for library in sample.libraries:
        logging.info("Analyzing library: %s" % (library.id)) 
        #
        # create library directory
        # 
        if not os.path.exists(library.output_dir):
            logging.info("Creating directory: %s" % (library.output_dir))
            os.makedirs(library.output_dir)
        deps.extend(run_library(library, genome, server, pipeline, num_processors, 
                                submit_job_func))
    return deps

def run_analysis(analysis_file, config_file, server_name,
                 num_processors, keep_tmp, rm_fastq):
    """
    rm_fastq: (True/False) delete fastq files after run
    keep_tmp: (True/False) delete temporary files after run
    """
    #
    # read configuration files
    #
    analysis = AnalysisConfig.from_xml(analysis_file)
    pipeline = PipelineConfig.from_xml(config_file)
    #
    # validate configuration files
    #
    if not analysis.is_valid():
        logging.error("Analysis config not valid")
        return config.JOB_ERROR
    if not pipeline.is_valid(server_name, analysis.species):
        logging.error("Pipeline config not valid")
        return config.JOB_ERROR
    # setup server
    server = pipeline.servers[server_name]
    if server.pbs:
        submit_job_func = submit_job_pbs
    else:
        submit_job_func = submit_job_nopbs
    # setup genome references
    genome = pipeline.species[analysis.species]
    #
    # attach analysis to pipeline output directory
    #
    analysis.attach_to_results(server.output_dir)    
    #
    # process each sample in analysis
    #
    deps = []
    for sample in analysis.samples:
        logging.info("Analyzing sample: %s" % (sample.id))        
        #
        # create output dir
        #
        if not os.path.exists(sample.output_dir):
            logging.info("Creating sample directory: %s" % (sample.output_dir))
            os.makedirs(sample.output_dir)
        #
        # process sample
        #        
        sample_deps = run_sample(sample, genome, server, pipeline, 
                                 num_processors, submit_job_func)
        #
        # cleanup temporary files
        #
        if not keep_tmp:
            msg = "Cleaning up temporary files"
            logging.info(msg)
            args = [sys.executable, os.path.join(_pipeline_dir, "cleanup_intermediate_files.py")]
            if rm_fastq:
                args.append("--rm-fastq")
            args.extend([os.path.abspath(analysis_file), server.output_dir, sample.id])
            job_id = submit_job_func("rm_%s" % (sample.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=sample.output_dir,
                                     walltime="1:00:00",
                                     deps=sample_deps)
            sample_deps = [job_id]
        #
        # write file indicating job is complete
        #
        msg = "Notifying user that job is complete"
        if os.path.exists(sample.job_complete_file) and (len(sample_deps) == 0):
            logging.info("[SKIPPED]: %s" % msg)
        else:
            logging.info(msg)
            args = [sys.executable, os.path.join(_pipeline_dir, "notify_complete.py"),
                    sample.job_complete_file]
            job_id = submit_job_func("done_%s" % (sample.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=sample.output_dir,
                                     walltime="1:00:00",
                                     email="ae",
                                     deps=sample_deps)
            sample_deps = [job_id]
        deps.extend(sample_deps)        
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--keep-tmp", action="store_true", dest="keep_tmp", default=False)
    parser.add_argument("--rm-fastq", action="store_true", dest="rm_fastq", default=False)
    parser.add_argument("analysis_file")
    parser.add_argument("config_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    return run_analysis(args.analysis_file, args.config_file, 
                        args.server_name, 
                        num_processors=args.num_processors,
                        keep_tmp=args.keep_tmp,
                        rm_fastq=args.rm_fastq)
    
if __name__ == '__main__': 
    sys.exit(main())