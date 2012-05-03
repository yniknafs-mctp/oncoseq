'''
Created on Mar 15, 2012

@author: oabalbin
'''

import argparse
import logging
import subprocess
import sys
import os

from oncoseq.lib import config
from oncoseq.lib.config import AnalysisConfig, PipelineConfig
from oncoseq.lib.cluster import submit_job_pbs, submit_job_nopbs
from oncoseq.lib.base import up_to_date

# setup pipeline script files
import oncoseq.lib
import oncoseq.exome.pipeline
import oncoseq.rnaseq.pipeline
_lib_dir=oncoseq.lib.__path__[0]
_pipeline_dir = oncoseq.exome.pipeline.__path__[0] 
_pipeline_rnaseq_dir = oncoseq.rnaseq.pipeline.__path__[0] 

#TODO: Adjust pmem= int(round(float(server.node_mem)/2, 0)), according to the input file size.
# This way I am asking 22 GB of memory which could be overkilling.
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
            pyscript = os.path.join(_lib_dir, "copy_uncompress_fastq.py")    
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
    qual_fastq_deps=[]
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
        qual_fastq_deps.extend([job_id])


    # TODO:
    # align against contaminants and filter out contaminant reads 
    #
    

    #
    # Align sequences using BWA
    #

    aln_job_ids=[]
    aln_sai_files=[]
    bwa_dep=[]
    skip_run = True

    for readnum in xrange(len(lane.copied_fastq_files)):    
        skip_run = skip_run and up_to_date(lane.exome_sam_aln,lane.copied_fastq_files[readnum]) \
        and up_to_date(lane.fastqc_report_files[readnum], lane.copied_fastq_files[readnum])
    msg = "Running BWA Alignment of exome sequences"
    if skip_run:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info("%s: Running BWA Alignment" % (lane.id))
        for readnum in xrange(len(lane.copied_fastq_files)):
            fastq_file = lane.copied_fastq_files[readnum]
            print fastq_file, readnum   
            logging.info('%s: calling BWA aln for mate %d fastq file %s' % (lane.id, readnum, fastq_file))
            
            aln_read = os.path.join(lane.output_dir, 'mate%d.sai' % (readnum))
            if up_to_date(aln_read, fastq_file):
                logging.info("[SKIPPED] aln file for mate%d.sai is up to date" % (readnum))
                aln_sai_files.append(aln_read)
            else:
                args = [pipeline.bwa_bin, "aln","-k",pipeline.bwa_config.nmismatch,"-t",pipeline.bwa_config.bwa_cores,
                        "-q",pipeline.bwa_config.bwa_qual_trim,"-f",aln_read,
                        os.path.join(server.references_dir, genome.get_path("genome_bwa_index")), 
                        fastq_file]
                
                # It runs with multi thread. Usually I run it with 4 or 6 processors.
                log_file = os.path.join(log_dir, "bwa_aln_reads%d.log" % (readnum+1))
                logging.debug("\targs: %s" % (' '.join(map(str, args))))
                job_id = submit_job_func("bwaA%d_%s" % (readnum+1, lane.id), args,
                                         num_processors=num_processors,
                                         node_processors=server.node_processors,
                                         node_memory=server.node_mem,
                                         pbs_script_lines=server.pbs_script_lines,
                                         working_dir=lane.output_dir,
                                         walltime="24:00:00",
                                         deps=qual_fastq_deps,
                                         stderr_filename=log_file)
                
                aln_sai_files.append(aln_read)
                aln_job_ids.extend([job_id])
                
        read_group="\'@RG\\tID:%s\\tPL:%s\\tLB:%s\\tSM:%s\'"%(lane.id, lane.platform, lane.library_id, lane.id)
        if len(aln_sai_files) == 1:
            
            args = [pipeline.bwa_bin, "samse","-f",lane.exome_sam_aln,
                    "-r",read_group,
                    os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                    aln_sai_files[0],lane.copied_fastq_files[0]]
        if len(aln_sai_files) == 2:
            args = [pipeline.bwa_bin, "sampe","-f",lane.exome_sam_aln,
                    "-r",read_group,
                    os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                    aln_sai_files[0],aln_sai_files[1],lane.copied_fastq_files[0],lane.copied_fastq_files[1]]
        
        # ask for 8 GB of memory by processor. It runs in 1 processor
        log_file = os.path.join(log_dir, "bwa_mapping_%s.log" % (lane.id))
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("bwaM_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pmem= int(round(float(server.node_mem)/2, 0)),
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="60:00:00",
                                 deps=aln_job_ids,
                                 stderr_filename=log_file)
        bwa_dep.extend([job_id])

                
    #
    # Convert SAM to BAM and limit the output to unique reads
    #
    sam_deps=[]
    skip_run = True
    skip_run = skip_run and up_to_date(lane.exome_bam_aln,lane.exome_sam_aln) \
        and up_to_date(lane.fastqc_report_files[readnum], lane.copied_fastq_files[readnum])
    msg = "Running SAM to BAM conversion of aligned sequences"
    if skip_run:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info("%s of lane %s" % (msg,lane.id))
        # lane.exome_bam_tmp
        #tmp_bam=open(lane.exome_bam_tmp,'w')
        #tmp_bam.close()        
        args = [pipeline.samtools_bin, "view", "-bSq", pipeline.bwa_config.mapping_qual,
                lane.exome_sam_aln, "-o",lane.exome_bam_tmp]
        log_file = os.path.join(log_dir, "sam2bam_%s.log" % (lane.id))
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        
        job_id = submit_job_func("sam2bam_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pmem= int(round(float(server.node_mem)/2, 0)),
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="24:00:00",
                                 deps=bwa_dep,
                                 stderr_filename=log_file)
        
        sam_deps.extend([job_id])

    #
    # Sort BAM file
    #
    
    sam_sort_deps=[]    
    msg = "Sorting BAM file of mapped reads"
    if up_to_date(lane.exome_bam_sorted, lane.exome_bam_aln):    
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "SortSam.jar"),
                "INPUT=%s" % (lane.exome_bam_tmp),
                "OUTPUT=%s" % (lane.exome_bam_sorted),
                "SO=coordinate",
                "CREATE_INDEX=true",
                "TMP_DIR=%s" % tmp_dir,
                "VALIDATION_STRINGENCY=SILENT"] 
        # TODO: Fix the following problem
        # "VALIDATION_STRINGENCY=SILENT" is necessary because "Exception in thread "main" 
        #net.sf.samtools.SAMFormatException: SAM validation error: ERROR: Record 6934, Read name 
        #PATHBIO-SOLEXA2_30LEJAAXX:7:1:82:1865, MAPQ should be 0 for unmapped read"
        
        log_file = os.path.join(log_dir, "sort_aligned_reads.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("sortexom_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pmem= int(round(float(server.node_mem)/2, 0)),
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="12:00:00",
                                 deps=sam_deps,
                                 stderr_filename=log_file)
        
        sam_sort_deps.extend([job_id])
        
        args=["rm",lane.exome_bam_tmp]
        msg = "Removing the temporary BAM file"
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("rm_tmp_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime="1:00:00",
                                 deps=sam_sort_deps,
                                 stderr_filename=log_file)
    

    #
    # Make alignment report
    #
        
    return sam_sort_deps



def run_library(library, genome, server, pipeline, num_processors,
                submit_job_func):
    
    
    all_lane_deps = []    
    library_deps=[]
    for lane in library.lanes:
        logging.info("Analyzing lane: %s" % (lane.id))
        
        #
        
        #
        # create lane directory
        # 
        if not os.path.exists(lane.output_dir):
            logging.info("Creating directory: %s" % (lane.output_dir))
            os.makedirs(lane.output_dir)
        #
        # process lane
        #
        thlane_deps=run_lane(lane, genome, server, pipeline, num_processors,submit_job_func)
        all_lane_deps.extend(thlane_deps)
    
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
    # merge lane alignment results
    #
    # If there is more than one lane for library, merge lane
    # files. 
    merge_bam_deps = []
    if len(library.lanes)!=1:
        msg = "Merging lane BAM files"
        if all(up_to_date(library.merged_bam_efile, lane.exome_bam_sorted) for lane in library.lanes):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            args = ["java", "-jar", 
                    os.path.join(pipeline.picard_dir, "MergeSamFiles.jar")]
            for lane in library.lanes:
                args.append("INPUT=%s" % (lane.tophat_bam_file))
            args.extend(["OUTPUT=%s" % (library.merged_bam_efile),
                         "SO=coordinate",
                         "MERGE_SEQUENCE_DICTIONARIES=true",
                         "USE_THREADING=true",
                         "CREATE_INDEX=true",
                         "TMP_DIR=%s" % tmp_dir,
                         "VALIDATION_STRINGENCY=SILENT"])
            log_file = os.path.join(log_dir, "picard_merge_bam_files.log")
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("merge_%s" % (library.id), args,
                                     num_processors=min(1, num_processors),
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pmem=int(round(float(server.node_mem)/2, 0)),
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=library.output_dir,
                                     walltime="24:00:00",
                                     deps=all_lane_deps,
                                     stderr_filename=log_file)
            merge_bam_deps = [job_id]
    else:
        msg = "[SKIPPED] Merging lane BAM files because only one lane for library"
        logging.info(msg)
        log_file = os.path.join(log_dir, "picard_merge_bam_files.log")
        for lane in library.lanes:
            lane_bam_index=lane.exome_bam_sorted.replace('.bam','.bai')
            lib_bam_index=library.merged_bam_efile.replace('.bam','.bai')
            source_files=[lane.exome_bam_sorted,lane_bam_index]
            target_files=[library.merged_bam_efile,lib_bam_index]
            for source, target in zip(source_files,target_files):
                args=["ln","-s",source, target]
                logging.debug("\targs: %s" % (' '.join(map(str, args))))
                job_id = submit_job_func("ln_%s" % (library.id), args,
                                         num_processors=min(1, num_processors),
                                         node_processors=server.node_processors,
                                         node_memory=server.node_mem,
                                         pbs_script_lines=server.pbs_script_lines,
                                         working_dir=library.output_dir,
                                         walltime="1:00:00",
                                         deps=all_lane_deps,
                                         stderr_filename=log_file)
    
                merge_bam_deps = [job_id] 
    
    #
    # Clean the BAM file before calling SNVs
    #
    msg = "Cleaning BAM Files"
    cleaning_dep = []
    if up_to_date(library.merged_cleaned_bam_efile, library.merged_bam_efile):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_lib_dir, "bam_cleaner.py"),
                "--picard-dir",pipeline.picard_dir,
                "--tmp-dir",tmp_dir,
                library.merged_bam_efile,
                library.merged_cleaned_bam_efile]
        
        log_file = os.path.join(log_dir, "bam_cleaning.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("bamcln_%s" % (library.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pmem=int(round(float(server.node_mem)/2, 0)),
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime="60:00:00",
                                 deps=merge_bam_deps,
                                 stderr_filename=log_file)
        cleaning_dep = [job_id]
    
    library_deps.extend(cleaning_dep)
    
    return library_deps


def run_sample(sample, genome, server, pipeline, num_processors,
               submit_job_func):
    sample_deps = []
    library_deps = []
    for library in sample.libraries:
        logging.info("Analyzing library: %s" % (library.id)) 
        #
        # create library directory
        # 
        if not os.path.exists(library.output_dir):
            logging.info("Creating directory: %s" % (library.output_dir))
            os.makedirs(library.output_dir)
        library_deps.extend(run_library(library, genome, server, pipeline, num_processors, 
                                submit_job_func))
        
    # TODO Merging at the sample level if more than one exome library by sample
    
        # 
    # create directories for tmp and log files
    #
    tmp_dir = os.path.join(sample.output_dir, "tmp") 
    if not os.path.exists(tmp_dir):
        logging.info("Creating directory: %s" % (tmp_dir))
        os.makedirs(tmp_dir)
    log_dir = os.path.join(sample.output_dir, "log")
    if not os.path.exists(log_dir):
        logging.info("Creating directory: %s" % (log_dir))
        os.makedirs(log_dir)
    
    #
    # merge lane alignment results
    #
    # If there is more than one lane for library, merge lane
    # files. 
    merge_bam_deps = []
    if len(sample.libraries)!=1:
        msg = "Merging lane BAM files"
        if all(up_to_date(sample.merged_bam_efile, library.merged_cleaned_bam_efile) for library in sample.libraries):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            args = ["java", "-jar", 
                    os.path.join(pipeline.picard_dir, "MergeSamFiles.jar")]
            for library in sample.libraries:
                args.append("INPUT=%s" % (library.merged_cleaned_bam_efile))
            args.extend(["OUTPUT=%s" % (sample.merged_bam_efile),
                         "SO=coordinate",
                         "MERGE_SEQUENCE_DICTIONARIES=true",
                         "USE_THREADING=true",
                         "CREATE_INDEX=true",
                         "TMP_DIR=%s" % tmp_dir,
                         "VALIDATION_STRINGENCY=SILENT"])
            log_file = os.path.join(log_dir, "picard_merge_bam_files.log")
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("merge_%s" % (library.id), args,
                                     num_processors=min(1, num_processors),
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pmem=int(round(float(server.node_mem)/2, 0)),
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=library.output_dir,
                                     walltime="24:00:00",
                                     deps=library_deps,
                                     stderr_filename=log_file)
            merge_bam_deps = [job_id]
            #
            # Clean the BAM file before calling SNVs
            #
            msg = "Cleaning BAM Files"
            cleaning_dep = []
            if up_to_date(sample.merged_cleaned_bam_efile, sample.merged_bam_efile):
                logging.info("[SKIPPED] %s" % msg)
            else:
                logging.info(msg)
                args = [sys.executable, os.path.join(_lib_dir, "bam_cleaner.py"),
                        "--picard-dir",pipeline.picard_dir,
                        "--tmp-dir",tmp_dir,
                        sample.merged_bam_efile,
                        sample.merged_cleaned_bam_efile]
                
                log_file = os.path.join(log_dir, "bam_cleaning.log")
                logging.debug("\targs: %s" % (' '.join(map(str, args))))
                job_id = submit_job_func("bamcln_%s" % (library.id), args,
                                         num_processors=1,
                                         node_processors=server.node_processors,
                                         node_memory=server.node_mem,
                                         pmem=int(round(float(server.node_mem)/2, 0)),
                                         pbs_script_lines=server.pbs_script_lines,
                                         working_dir=library.output_dir,
                                         walltime="60:00:00",
                                         deps=merge_bam_deps,
                                         stderr_filename=log_file)
                cleaning_dep = [job_id]

    else:
        msg = "[SKIPPED] Merging lane BAM files because only one library for sample"
        logging.info(msg)
        log_file = os.path.join(log_dir, "picard_merge_bam_files.log")        
        for library in sample.libraries:
            lib_bam_index=library.merged_cleaned_bam_efile.replace('.bam','.bai')
            sample_bam_index=sample.merged_cleaned_bam_efile.replace('.bam','.bai')            
            source_files=[library.merged_cleaned_bam_efile,lib_bam_index]
            target_files=[sample.merged_cleaned_bam_efile,sample_bam_index]
            # create symbolic links        
            for source, target in zip(source_files,target_files):
                args=["ln","-s",source, target]
                logging.debug("\targs: %s" % (' '.join(map(str, args))))
                
                job_id = submit_job_func("ln_%s" % (library.id), args,
                                         num_processors=min(1, num_processors),
                                         node_processors=server.node_processors,
                                         node_memory=server.node_mem,
                                         pbs_script_lines=server.pbs_script_lines,
                                         working_dir=library.output_dir,
                                         walltime="1:00:00",
                                         deps=library_deps,
                                         stderr_filename=log_file)

                cleaning_dep = [job_id]
    

    
    # 
    # Determine homozygous status for cosmic postions in actionable genes
    #
    
    msg = "Determining homozygous status for cosmic postions in actionable genes"
    homo_dep = [] 
    
    if up_to_date(sample.cosmic_qual_vcf_file, sample.merged_cleaned_bam_efile):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_pipeline_dir, "cosmic_positions_homozygosity.py"),
                os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                sample.merged_cleaned_bam_efile,
                "None", # Any samtools_bcf_file not used
                sample.cosmic_qual_vcf_file,
                os.path.join(server.references_dir, genome.get_path("cosmic_positions"))]
        log_file = os.path.join(log_dir, "cosmic_postions_homozygosity.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("cosmicqual_%s" % (sample.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime="24:00:00",
                                 deps=cleaning_dep,
                                 stderr_filename=log_file)
        homo_dep = [job_id]    
    
    
    # TODO
    # Compute Exon Coverage
    #
    '''
    msg = "Compute exome coverage"
    cov_pro_dep = [] 
    if up_to_date(sample.exome_coverage_file, sample.merged_cleaned_bam_efile):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_pipeline_dir, "target_coverage.py"),
                os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                sample.merged_cleaned_bam_efile,
                "", # Any samtools_bcf_file not used
                sample.exome_coverage_file,
                os.path.join(server.references_dir, genome.get_path("exome_bed_file"))]
        log_file = os.path.join(log_dir, "exome_coverage.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("probe_cov_%s" % (sample.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime="60:00:00",
                                 deps=cleaning_dep,
                                 stderr_filename=log_file)
        cov_pro_dep = [job_id]    
    '''
    # TODO
    # Compute Probe Coverage
    #
    
    msg = "Compute probe coverage"
    cov_pro_dep = [] 
    if up_to_date(sample.probe_summary_file, sample.merged_cleaned_bam_efile):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_pipeline_dir, "target_coverage.py"),
                sample.merged_cleaned_bam_efile,
                os.path.join(server.references_dir, genome.get_path("capture_roche")),
                sample.probe_coverage_file,
                sample.probe_summary_file,
                pipeline.vscan_config.min_avgbase_quality]
        
        log_file = os.path.join(log_dir, "capture_coverage.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("probe_cov_%s" % (sample.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime="24:00:00",
                                 deps=cleaning_dep,
                                 stderr_filename=log_file)
        cov_pro_dep = [job_id]    
    
    #
    # generate Exome coverage bedgraph  
    #
    msg = "Generating coverage bedgraph file"
    bedgraph_deps = []
    if (up_to_date(sample.coverage_bedgraph_file, sample.merged_cleaned_bam_efile)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [os.path.join(pipeline.bedtools_dir, "genomeCoverageBed"),
                "-bg", "-split", "-g", (os.path.join(server.references_dir, genome.get_path("chrom_sizes"))),
                "-ibam", sample.merged_cleaned_bam_efile]
        log_file = os.path.join(log_dir, "genomeCoverageBed.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("bedtools_%s" % (sample.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime="20:00:00",
                                 deps=cleaning_dep,
                                 stdout_filename=sample.coverage_bedgraph_file,
                                 stderr_filename=log_file)
        bedgraph_deps = [job_id]
    #
    # convert bedgraph to bigwig coverage file
    #
    msg = "Create bigWig file to display coverage"
    bigwig_deps = []
    if (up_to_date(sample.coverage_bigwig_file, sample.coverage_bedgraph_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [os.path.join(pipeline.ucsc_dir, "bedGraphToBigWig"),
                sample.coverage_bedgraph_file,
                os.path.join(server.references_dir, genome.get_path("chrom_sizes")),
                sample.coverage_bigwig_file]
        log_file = os.path.join(log_dir, "bedGraphToBigWig.log")
        job_id = submit_job_func("bigwig_%s" % (sample.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime="20:00:00",
                                 deps=bedgraph_deps,
                                 stderr_filename=log_file)
        bigwig_deps = [job_id]

    

    deps=cleaning_dep + homo_dep + cov_pro_dep + bigwig_deps #+ cov_pro_dep #+ cov_exon_dep
    sample_deps.extend(deps)
            
    return sample_deps


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
    # 
    # create directories for tmp and log files
    #    
    patient_deps = []
    for patient, patient_samples in analysis.patients.iteritems():
        
        for sample in patient_samples:
            if sample.protocol != config.EXOME:
                logging.info("IN EXOME protocol SKIPPING sample: %s, because protocol is %s" % (sample.id,sample.protocol))
                continue
            # TODO: This could be changed by class A and class B
            if sample.cancer_progression == "benign":
                benign_sample = sample.merged_cleaned_bam_efile
                benign_cov = sample.probe_summary_file
                    
            else:
                tumor_sample = sample.merged_cleaned_bam_efile
                tumor_cov = sample.probe_summary_file
                
                
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
        
        # Finish work at sample level for this patient
        patient_deps.extend(sample_deps)
        
        tmp_dir = os.path.join(patient.output_dir, "tmp") 
        if not os.path.exists(tmp_dir):
            logging.info("Creating directory: %s" % (tmp_dir))
            os.makedirs(tmp_dir)
        log_dir = os.path.join(patient.output_dir, "log")
        if not os.path.exists(log_dir):
            logging.info("Creating directory: %s" % (log_dir))
            os.makedirs(log_dir)
        
        # TODO: Analysis at the patient level:
        # Calling of somatic mutations
        # Calling of Somatic CNVs
        #       Calling of LOH
        #       Calling of Germline Mutations
        # Determine confidence for particular cosmic positions
        # Compute coverage over the probe regions
        
        
        #
        # call somatic snps
        #
        msg = "Calling SOMATIC SNVs with samtools"
        samtools_snv_deps = [] 
        
        if up_to_date(patient.samtools_vcf_file, sample.merged_cleaned_bam_efile):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            args = [sys.executable, os.path.join(_pipeline_dir, "snps_bySamtools_Somatic.py"),
                    os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                    benign_sample,
                    tumor_sample,
                    patient.samtools_bcf_file,
                    patient.samtools_vcf_file]
            log_file = os.path.join(log_dir, "samtools_snp_calling.log")
            log_stderr_file = os.path.join(log_dir, "samtools_snp_calling_stderr.log")
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("samsnv_%s" % (patient.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pmem= int(round(float(server.node_mem)/2, 0)),
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=patient.output_dir,
                                     walltime="60:00:00",
                                     deps=patient_deps,
                                     stderr_filename=log_file)
            samtools_snv_deps = [job_id]    
        # 
        # Call SNVs using varscan
        #
                
        msg = "Calling SOMATIC SNVs with VarScan"
        varscan_deps = []
        if up_to_date(patient.varscan_snv_file, sample.merged_cleaned_bam_efile):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            args = [sys.executable, os.path.join(_pipeline_dir, "snps_byVarScan_Somatic.py"),
                    "--varscan-dir", pipeline.varscan_dir,
                    os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                    benign_sample,
                    tumor_sample,
                    patient.varscan_snv_file,
                    patient.varscan_indel_file]
            log_stdout_file = os.path.join(log_dir, "varscan_snp_calling_stdout.log")
            log_stderr_file = os.path.join(log_dir, "varscan_snp_calling_stderr.log")
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("varscan_%s" % (patient.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pmem= int(round(float(server.node_mem)/2, 0)),
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=patient.output_dir,
                                     walltime="60:00:00",
                                     deps=patient_deps,
                                     stdout_filename=log_stdout_file,
                                     stderr_filename=log_stderr_file)
            varscan_deps = [job_id]

        call_deps=samtools_snv_deps + varscan_deps
        
        
        #
        # Running CNV analysis with Exome CNV
        #
        # Determining the Read length for the tumor sample
        # TODO: This could be modified into a particular small function
        '''
        args = [pipeline.samtools_bin,"view",tumor_sample,"|", "head", "-n", "1", "|", "cut", "-f", "10"]                
        args=",".join(args).replace(',',' ')
        p1 =subprocess.Popen(args,stdout=subprocess.PIPE,shell=True)
        t = p1.communicate()[0]
        tumor_read_length=len(t.strip('\n'))
        '''
        msg = "Calling CNVs with Exome CNVs"
        cnv_deps = []
        if up_to_date(patient.exome_cnv_file, sample.merged_cleaned_bam_efile):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            args = [sys.executable, os.path.join(_pipeline_dir, "exomeCNV.py"),
                    pipeline.samtools_bin, 
                    pipeline.rscript_bin,
                    benign_cov,
                    tumor_cov,
                    patient.exome_loh_file,
                    patient.exome_cnv_file,
                    patient.exome_cnv_plot,
                    tumor_sample]#,"TRUE"]
            
            log_stderr_file = os.path.join(log_dir, "exome_cnv_calling_stderr.log")
            log_stdout_file = os.path.join(log_dir, "exome_cnv_calling_stdout.log")
            logging.debug("\targs: %s" % (' '.join(map(str, args))))
            job_id = submit_job_func("cnv_%s" % (patient.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=patient.output_dir,
                                     walltime="24:00:00",
                                     deps=patient_deps,
                                     stdout_filename=log_stdout_file,
                                     stderr_filename=log_stderr_file)
            cnv_deps = [job_id]

        call_deps.extend(cnv_deps)
        
        #
        # write file indicating patient job is complete
        #
        msg = "Notifying user that job is complete"
        if os.path.exists(patient.job_complete_file) and (len(sample_deps) == 0):
            logging.info("[SKIPPED]: %s" % msg)
        else:
            logging.info(msg)
            args = [sys.executable, os.path.join(_lib_dir, "notify_complete.py"),
                    patient.job_complete_file]
            job_id = submit_job_func("done_%s" % (patient.id), args,
                                     num_processors=1,
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=sample.output_dir,
                                     walltime="1:00:00",
                                     email="ae",
                                     deps=call_deps) # TODO: patient dependency
            
        call_deps.extend([job_id])
        
            
                  
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
