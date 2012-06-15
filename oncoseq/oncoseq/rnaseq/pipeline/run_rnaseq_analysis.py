'''
Created on Aug 7, 2011

@author: mkiyer
'''
import sys
import os
import logging
import xml.etree.cElementTree as etree

from oncoseq.lib import config
from oncoseq.lib.base import up_to_date, indent_xml
from oncoseq.lib.seqdb import SAMPLE_TYPE_RNASEQ, SAMPLE_TYPE_CAPTURE_RNASEQ

# setup pipeline script files
import oncoseq.pipeline
import oncoseq.rnaseq.pipeline
_oncoseq_pipeline_dir = oncoseq.pipeline.__path__[0]
_rnaseq_pipeline_dir = oncoseq.rnaseq.pipeline.__path__[0] 
num_processes=3 # tmp move to the config file.

def run_lane(lane, genome, server, pipeline, num_processors,
             submit_job_func):
    #
    # create lane directory
    # 
    if not os.path.exists(lane.output_dir):
        logging.info("Creating directory: %s" % (lane.output_dir))
        os.makedirs(lane.output_dir)
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
    # run fastqc program
    #
    fastqc_deps = []
    skip_run = True
    for readnum in xrange(len(lane.fastq_files)):    
        skip_run = skip_run and (up_to_date(lane.fastqc_data_files[readnum], lane.fastq_files[readnum]) and
                                 up_to_date(lane.fastqc_report_files[readnum], lane.fastq_files[readnum]))
    msg = "Running FASTQC quality assessment"
    if skip_run:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info("%s" % msg)       
        num_threads = min(num_processors, len(lane.fastq_files))
        args = [pipeline.fastqc_bin, 
                "--threads", num_threads,
                "-o", lane.output_dir]
        args.extend(lane.fastq_files)        
        log_file = os.path.join(log_dir, "fastqc.log")
        job_id = submit_job_func("fqc_%s" % (lane.id), args,
                                 num_processors=num_threads,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="10:00:00",
                                 stderr_filename=log_file)
        fastqc_deps.append(job_id)
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
            pyscript = os.path.join(_oncoseq_pipeline_dir, "copy_uncompress_fastq.py")    
            args = [sys.executable, pyscript,
                    "--quals", lane.quality_scores,
                    "--fastqc-data-files", ",".join(lane.fastqc_data_files),
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
                                     deps=fastqc_deps,
                                     stderr_filename=log_file)
            copy_fastq_deps.append(job_id)
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
        args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "filter_abundant_sequences.py"),
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
        args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "filter_unmapped_pairs.py"),
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
                os.path.join(_rnaseq_pipeline_dir, "estimate_fragment_size_distribution.py"),
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
        args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "run_tophat.py"),
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
                     os.path.join(server.references_dir, genome.get_path("genome_bowtie2_index")),
                     lane.frag_size_dist_file])
        args.extend(lane.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        log_file = os.path.join(log_dir, "tophat.log")
        # allocate 16gb to run tophat
        job_id = submit_job_func("tophat_%s" % (lane.id), args,
                                 num_processors=num_processors,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 mem=16000,
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
    bigwig_deps = []
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
        bigwig_deps = [job_id]
    return tophat_deps, tophat_deps + bigwig_deps


def run_library(library, genome, server, pipeline, num_processors,
                submit_job_func):
    # create library directory
    if not os.path.exists(library.output_dir):
        logging.info("Creating directory: %s" % (library.output_dir))
        os.makedirs(library.output_dir)       
    merge_lane_deps = []
    all_lane_deps = []
    has_paired_end = False
    frag_size_dist_files = []
    for lane in library.lanes:
        logging.info("Analyzing lane: %s" % (lane.id))
        #
        # process lane
        #
        tophat_deps, lane_deps = run_lane(lane, genome, server, pipeline, num_processors, submit_job_func)
        merge_lane_deps.extend(tophat_deps)
        all_lane_deps.extend(lane_deps)
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
                os.path.join(_rnaseq_pipeline_dir, "merge_fragment_size_distributions.py"),
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
                                 deps=merge_lane_deps,
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
                                 deps=merge_lane_deps,
                                 stderr_filename=log_file)
        merge_bam_deps = [job_id]
    
    
    #
    # Clean the BAM file before calling SNVs
    #
    msg = "Cleaning BAM Files"
    cleaning_dep = []
    if up_to_date(library.merged_cleaned_bam_file, library.merged_bam_file):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_oncoseq_pipeline_dir, "bam_cleaner.py"),
                "--picard-dir",pipeline.picard_dir,
                "--tmp-dir",tmp_dir,
                library.merged_bam_file,
                library.merged_cleaned_bam_file]       
        log_file = os.path.join(log_dir, "bam_cleaning.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("bamcln_%s" % (library.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="60:00:00",
                                 deps=merge_bam_deps,
                                 stderr_filename=log_file)
        cleaning_dep = [job_id]

    #
    # call snps
    #
    msg = "Calling SNVs with samtools"
    samtools_snv_deps = [] 
    
    if up_to_date(library.samtools_vcf_file, library.merged_bam_file):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "call_snps.py"),
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
                                 deps=cleaning_dep,
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
        args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "call_snps_varscan.py"),
                "--varscan-dir", pipeline.varscan_dir,
                os.path.join(server.references_dir, genome.get_path("genome_fasta_file")),
                library.merged_bam_file,
                library.varscan_snv_file,
                library.varscan_indel_file]
        log_stdout_file = os.path.join(log_dir, "varscan_snp_calling_stdout.log")
        log_stderr_file = os.path.join(log_dir, "varscan_snp_calling_stderr.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("varscan_%s" % (library.id), args,
                                 num_processors=server.node_processors/num_processes,#1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 mmem= int(round(float(server.node_mem)/num_processes, 0)), #int(round(float(server.node_mem)/2, 0)),
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="60:00:00",
                                 deps=cleaning_dep,
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
        args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "run_cufflinks.py"),
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
    lib_deps = cufflinks_deps + samtools_snv_deps + varscan_deps
    return lib_deps + all_lane_deps

def run_sample(sample, genome, server, pipeline, num_processors,
               submit_job_func, keep_tmp):
    # create directory
    if not os.path.exists(sample.output_dir):
        logging.info("Creating directory: %s" % (sample.output_dir))
        os.makedirs(sample.output_dir)
    # output sample XML
    root = etree.Element("analysis")
    sample.to_xml(root)
    indent_xml(root)
    f = open(sample.xml_file, "w")
    print >>f, etree.tostring(root)
    f.close()
    # run libraries
    lib_deps = []
    for library in sample.libraries:
        logging.info("Analyzing library: %s" % (library.id)) 
        lib_deps.extend(run_library(library, genome, server, pipeline, num_processors, 
                                    submit_job_func))
    #
    # write file indicating job is complete
    #
    deps = lib_deps
    msg = "Notifying user that job is complete"
    if os.path.exists(sample.job_complete_file) and (len(lib_deps) == 0):
        logging.info("[SKIPPED]: %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_oncoseq_pipeline_dir, "notify_complete.py"),
                sample.job_complete_file]
        job_id = submit_job_func("done_%s" % (sample.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime="1:00:00",
                                 email="ae",
                                 deps=lib_deps)
        deps = [job_id]
    return deps

def run_sample_group(grp, genome, server, pipeline, num_processors,
                     submit_job_func, keep_tmp):
    #
    # process sample
    #
    sample_deps = []
    for sample_type in (SAMPLE_TYPE_RNASEQ, SAMPLE_TYPE_CAPTURE_RNASEQ):        
        sample = grp.samples[sample_type]
        if sample is None:
            continue
        deps = run_sample(sample, genome, server, pipeline, 
                          num_processors, submit_job_func, keep_tmp)
        # cleanup temporary files
        if keep_tmp:
            sample_deps.extend(deps)
        else:
            msg = "Cleaning up temporary files"
            logging.info(msg)
            args = [sys.executable, os.path.join(_rnaseq_pipeline_dir, "cleanup_intermediate_files.py")]
            args.extend([os.path.abspath(sample.xml_file), grp.output_dir, sample.id])
            job_id = submit_job_func("rm_%s" % (sample.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=sample.output_dir,
                                     walltime="1:00:00",
                                     deps=deps)
            sample_deps.append(job_id)
    #
    # write file indicating job is complete
    #
    deps = []
    msg = "Notifying user that job is complete"
    if os.path.exists(grp.job_complete_file) and (len(sample_deps) == 0):
        logging.info("[SKIPPED]: %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_oncoseq_pipeline_dir, "notify_complete.py"),
                grp.job_complete_file]
        job_id = submit_job_func("done_%s" % (grp.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=grp.output_dir,
                                 walltime="1:00:00",
                                 email="ae",
                                 deps=sample_deps)
        deps = [job_id]
    return deps
