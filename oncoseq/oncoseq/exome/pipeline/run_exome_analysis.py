'''
Created on Mar 15, 2012

@author: oabalbin
'''
import logging
import sys
import os

from oncoseq.lib import config
from oncoseq.lib.base import up_to_date
from oncoseq.lib.seqdb import SAMPLE_TYPE_EXOME_TUMOR, SAMPLE_TYPE_EXOME_NORMAL

# setup pipeline script files
import oncoseq.pipeline
import oncoseq.exome.pipeline
_oncoseq_pipeline_dir = oncoseq.pipeline.__path__[0]
_exome_pipeline_dir = oncoseq.exome.pipeline.__path__[0] 

num_processes=3
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
    # run fastqc program
    #
    qual_fastq_deps = []
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
                                 mem=config.FASTQC_JOB_MEM,
                                 walltime=config.FASTQC_JOB_WALLTIME,
                                 stderr_filename=log_file)
        qual_fastq_deps.extend([job_id])
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
                                     mem=config.COPY_UNCOMPRESS_JOB_MEM,
                                     walltime=config.COPY_UNCOMPRESS_JOB_WALLTIME,
                                     deps=qual_fastq_deps,
                                     stderr_filename=log_file)
            copy_fastq_deps.append(job_id)

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
                                         mem=config.BWA_ALIGNMENT_JOB_MEM*num_processors,
                                         walltime=config.BWA_ALIGNMENT_JOB_WALLTIME,
                                         deps=copy_fastq_deps,
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
                                 mem=config.BWA_MAPPING_JOB_MEM,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime=config.BWA_MAPPING_JOB_WALLTIME,#"60:00:00",
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
                                 num_processors=server.node_processors/num_processes, #1
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 mem=config.SAM2BAM_JOB_MEM,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime=config.SAM2BAM_JOB_WALLTIME,
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
                                 num_processors=server.node_processors/num_processes,#1
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 mem= config.SORT_READS_JOB_MEM,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 walltime=config.SORT_READS_JOB_WALLTIME,
                                 deps=sam_deps,
                                 stderr_filename=log_file)
        
        sam_sort_deps.extend([job_id])
        # Removing Intermediate files and original sequence reads.
        args=["rm",lane.exome_bam_tmp, lane.exome_sam_aln]       
        for readnum in xrange(len(lane.copied_fastq_files)):
            fastq_file = lane.copied_fastq_files[readnum]
            args.append(fastq_file)
            
        msg = "Removing the temporary BAM/SAM files and copied sequence reads"
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("rm_tmp_%s" % (lane.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=lane.output_dir,
                                 mem=config.CLEANUP_INTERMEDIATE_FILES_JOB_MEM,
                                 walltime=config.CLEANUP_INTERMEDIATE_FILES_JOB_WALLTIME,
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
                                     num_processors=1, 
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     mem=config.MERGE_SAM_FILES_JOB_MEM,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=library.output_dir,
                                     walltime=config.MERGE_SAM_FILES_JOB_WALLTIME,
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
                                         mem=config.SYMBOLIC_LINK_JOB_MEM,
                                         walltime=config.SYMBOLIC_LINK_JOB_WALLTIME,
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
        args = [sys.executable, os.path.join(_oncoseq_pipeline_dir, "bam_cleaner.py"),
                "--picard-dir",pipeline.picard_dir,
                "--tmp-dir",tmp_dir,
                library.merged_bam_efile,
                library.merged_cleaned_bam_efile]
        
        log_file = os.path.join(log_dir, "bam_cleaning.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("bamcln_%s" % (library.id), args,
                                 num_processors=1, #1
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 mem=config.BAM_CLEANING_JOB_MEM,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=library.output_dir,
                                 walltime=config.BAM_CLEANING_JOB_WALLTIME,
                                 deps=merge_bam_deps,
                                 stderr_filename=log_file)
        cleaning_dep = [job_id]
    
    library_deps.extend(cleaning_dep)
    
    return library_deps


def run_sample(sample, genome, server, pipeline, num_processors,
               submit_job_func):
    # create output dir
    if not os.path.exists(sample.output_dir):
        logging.info("Creating sample directory: %s" % (sample.output_dir))
        os.makedirs(sample.output_dir)
    # process libraries
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
        
    # TODO: Merging at the sample level if more than one exome library by sample

    
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
                                     num_processors=1, #1
                                     node_processors=server.node_processors,
                                     node_memory=server.node_mem,
                                     mem=config.MERGE_SAM_FILES_JOB_MEM,
                                     pbs_script_lines=server.pbs_script_lines,
                                     working_dir=library.output_dir,
                                     walltime=config.MERGE_SAM_FILES_JOB_WALLTIME,
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
                args = [sys.executable, os.path.join(_oncoseq_pipeline_dir, "bam_cleaner.py"),
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
                                         mem=config.BAM_CLEANING_JOB_MEM,#int(round(float(server.node_mem)/num_processes, 0)), #pmem
                                         pbs_script_lines=server.pbs_script_lines,
                                         working_dir=library.output_dir,
                                         walltime=config.BAM_CLEANING_JOB_WALLTIME,
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
                                         mem=config.SYMBOLIC_LINK_JOB_MEM,
                                         walltime=config.SYMBOLIC_LINK_JOB_WALLTIME,
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
        args = [sys.executable, os.path.join(_exome_pipeline_dir, "cosmic_positions_homozygosity.py"),
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
                                 mem=config.COSMIC_COVERAGE_JOB_MEM,
                                 walltime=config.COSMIC_COVERAGE_JOB_WALLTIME,
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
        args = [sys.executable, os.path.join(_exome_pipeline_dir, "target_coverage.py"),
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
        args = [sys.executable, os.path.join(_exome_pipeline_dir, "target_coverage.py"),
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
                                 mem=config.CAPTURE_COVERAGE_JOB_MEM,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=sample.output_dir,
                                 walltime=config.CAPTURE_COVERAGE_JOB_WALLTIME,
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
                                 mem=config.COVERAGE_JOB_MEM,
                                 walltime=config.COVERAGE_JOB_WALLTIME,
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
                                 walltime=config.BIGWIG_JOB_MEM,
                                 mem=config.BIGWIG_JOB_WALLTIME,
                                 deps=bedgraph_deps,
                                 stderr_filename=log_file)
        bigwig_deps = [job_id]

    sample_deps=cleaning_dep + homo_dep + cov_pro_dep + bigwig_deps #+ cov_pro_dep #+ cov_exon_dep
    return sample_deps

def run_sample_group(grp, genome, server, pipeline, num_processors,
                     submit_job_func, keep_tmp):
    # check exome samples
    tumor_sample = grp.samples[SAMPLE_TYPE_EXOME_TUMOR]
    benign_sample = grp.samples[SAMPLE_TYPE_EXOME_NORMAL]
    if (benign_sample is None) or (tumor_sample is None):
        logging.info("Skipping exome analysis: tumor and/or benign exome samples missing")
        return config.JOB_SUCCESS
    # process samples
    sample_deps = run_sample(benign_sample, genome, server, pipeline, 
                             num_processors, submit_job_func)
    sample_deps.extend(run_sample(tumor_sample, genome, server, pipeline, 
                                  num_processors, submit_job_func))
    # make temp directories
    tmp_dir = os.path.join(grp.output_dir, "tmp") 
    if not os.path.exists(tmp_dir):
        logging.info("Creating directory: %s" % (tmp_dir))
        os.makedirs(tmp_dir)
    log_dir = os.path.join(grp.output_dir, "log")
    if not os.path.exists(log_dir):
        logging.info("Creating directory: %s" % (log_dir))
        os.makedirs(log_dir)
    #
    # Call somatic variants with samtools
    #
    msg = "Calling somatic SNVs with samtools"
    samtools_snv_deps = []
    if (up_to_date(grp.samtools_vcf_file, benign_sample.merged_cleaned_bam_efile) and
        up_to_date(grp.samtools_vcf_file, tumor_sample.merged_cleaned_bam_efile)):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_exome_pipeline_dir, "snps_bySamtools_Somatic.py"),
                os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                benign_sample.merged_cleaned_bam_efile,
                tumor_sample.merged_cleaned_bam_efile,
                grp.samtools_bcf_file,
                grp.samtools_vcf_file]
        log_file = os.path.join(log_dir, "samtools_snp_calling.log")
        log_stderr_file = os.path.join(log_dir, "samtools_snp_calling_stderr.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("samsnv_%s" % (grp.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 mem=config.SAMTOOLS_VARIANT_JOB_MEM,#pmem
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=grp.output_dir,
                                 walltime=config.SAMTOOLS_VARIANT_JOB_WALLTIME,
                                 deps=sample_deps,
                                 stderr_filename=log_file)
        samtools_snv_deps = [job_id]    
    # 
    # Call somatic variants with varscan
    #
    msg = "Calling somatic SNVs with Varscan"
    varscan_deps = []
    if (up_to_date(grp.varscan_snv_file, benign_sample.merged_cleaned_bam_efile) and
        up_to_date(grp.varscan_snv_file, tumor_sample.merged_cleaned_bam_efile)):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_exome_pipeline_dir, "snps_byVarScan_Somatic.py"),
                "--varscan-dir", pipeline.varscan_dir,
                os.path.join(server.references_dir, genome.get_path("genome_bwa_index")),
                benign_sample.merged_cleaned_bam_efile,
                tumor_sample.merged_cleaned_bam_efile,
                grp.varscan_snv_file,
                grp.varscan_indel_file]
        log_stdout_file = os.path.join(log_dir, "varscan_snp_calling_stdout.log")
        log_stderr_file = os.path.join(log_dir, "varscan_snp_calling_stderr.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("varscan_%s" % (grp.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 mem= config.VARSCAN_VARIANT_JOB_MEM,#2 pmem
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=grp.output_dir,
                                 walltime=config.VARSCAN_VARIANT_JOB_WALLTIME,
                                 deps=sample_deps,
                                 stdout_filename=log_stdout_file,
                                 stderr_filename=log_stderr_file)
        varscan_deps = [job_id]
    #
    # Running CNV analysis with Exome CNV
    #
    msg = "Calling CNVs with ExomeCNV"
    cnv_deps = []
    if (up_to_date(grp.exome_cnv_file, benign_sample.merged_cleaned_bam_efile) and
        up_to_date(grp.exome_cnv_file, tumor_sample.merged_cleaned_bam_efile)):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_exome_pipeline_dir, "exomeCNV.py"),
                pipeline.samtools_bin, 
                pipeline.rscript_bin,
                benign_sample.probe_summary_file,
                tumor_sample.probe_summary_file,
                grp.exome_loh_file,
                grp.exome_cnv_file,
                grp.exome_cnv_plot,
                tumor_sample.merged_cleaned_bam_efile]
        log_stderr_file = os.path.join(log_dir, "exome_cnv_calling_stderr.log")
        log_stdout_file = os.path.join(log_dir, "exome_cnv_calling_stdout.log")
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        job_id = submit_job_func("cnv_%s" % (grp.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=grp.output_dir,
                                 mem=config.EXOME_CNV_JOB_MEM,
                                 walltime=config.EXOME_CNV_JOB_WALLTIME,
                                 deps=sample_deps,
                                 stdout_filename=log_stdout_file,
                                 stderr_filename=log_stderr_file)
        cnv_deps = [job_id]
    # all dependencies for sample group
    grp_deps = samtools_snv_deps + varscan_deps + cnv_deps
    #
    # write file indicating sample group jobs are complete
    #
    deps = grp_deps
    msg = "Notifying user that sample group jobs are complete"
    if os.path.exists(grp.dna_job_complete_file) and (len(grp_deps) == 0):
        logging.info("[SKIPPED]: %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_oncoseq_pipeline_dir, "notify_complete.py"),
                grp.dna_job_complete_file]
        job_id = submit_job_func("dnadone_%s" % (grp.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=grp.output_dir,
                                 mem=config.NOTIFY_COMPLETE_JOB_MEM,
                                 walltime=config.NOTIFY_COMPLETE_JOB_WALLTIME,
                                 email="ae",
                                 deps=grp_deps)
        deps = [job_id]
    return deps
