'''
Created on Nov 27, 2010

@author: oabalbin
@author: mkiyer

Dec 8, 2010 -- [mkiyer] merged with main alignment pipeline
'''
import logging
import os
import sys
import subprocess
#import argparse
import tempfile
from collections import defaultdict, deque

from base import JOB_SUCCESS, JOB_ERROR, up_to_date
from job_runner import qsub_cac, qsub_loc, run_local
from config import ExomePipelineConfig, ExomeLaneConfig

# define base paths to job scripts
_jobs_dir = os.path.abspath(os.path.dirname(__file__))

def reads_alignment(ref_genome, reads_file_name, nmismatch, num_cores, outfile_name, path_to_bwa):
    '''
    time /exds/users/oabalbin/sw/bwa/bwa-0.5.8a/bwa aln -k 2 -t 3 
    /exds/projects/alignment_indexes/bwa/hg19_illumina/hg19_illumina.fa s_4_2_sequence.txt > 
    s_4_2_sequence.txt.sai
    '''
    quality_trim=str('20')
    bwa_command = os.path.join(path_to_bwa, 'bwa')
    
    # bwa.0.5.10 includes the option -I for illumina platform and -q could be included.
    # However, I prefered to convert the bases at the beginning of the pipeline so everything gets run with the same
    # files.
    
    #args=[bwa_command,'aln', '-k', str(nmismatch), '-t',str(num_cores), '-f', outfile_name, ref_genome, reads_file_name] #'-q', quality_trim
    #args=[bwa_command,'aln', '-k', str(nmismatch), '-t', str(num_cores), '-q', quality_trim, '-I', '-f', outfile_name, ref_genome, reads_file_name]
    args=[bwa_command,'aln', '-k', str(nmismatch), '-t', str(num_cores), '-q', quality_trim, '-f', outfile_name, ref_genome, reads_file_name]
    
    return args

def single_read_alignment_sam(ref_genome, align_reads_pair1, raw_reads1, outfile_name, path_to_bwa,
                              rg_id, rg_sample, rg_lib, rg_platform):
    '''
    time /exds/users/oabalbin/sw/bwa/bwa-0.5.8a/bwa sampe 
    /exds/projects/alignment_indexes/bwa/hg19_illumina/hg19_illumina.fa 
    ./s_4_1_sequence.txt.sai ./s_4_2_sequence.txt.sai ./s_4_1_sequence.txt 
    ./s_4_2_sequence.txt > s_4_12_sequence.aln.sam
    
    # Using the path for the gatk well formed bam files with read group, 
    and platform information
    time /exds/users/oabalbin/sw/bwa/bwa-0.5.8a/bwa sampe 
    -i MARK -m DEPRISTO -l LIBRARY -p ILLUMINA 
    /exds/projects/alignment_indexes/gatk/hg19/bwa/hg19.fa 
    ./s_3_1_sequence.hg19.sai ./s_3_2_sequence.hg19.sai 
    ./s_3_1_sequence.txt ./s_3_2_sequence.txt > s_3_12_sequence.hg19.aln2.sam

    bwa samse [-n max_occ] [-f out.sam] <prefix> <in.sai> <in.fq>

    BWA parameters for generating read group information:
         -i       read group identifier (ID)
         -m       read group sample (SM), required if ID is given
         -l       read group library (LB)
         -p       read group platform (PL)
    '''
    bwa_command = os.path.join(path_to_bwa, 'bwa')
    args = [bwa_command, 'samse', '-f', outfile_name, 
            '-i', rg_id, '-m', rg_sample, '-l', rg_lib, '-p', rg_platform,
            ref_genome, align_reads_pair1, raw_reads1]
    return args

def paired_read_alignment_sam(ref_genome, align_reads_pair1, align_reads_pair2, 
                              raw_reads1, raw_reads2, outfile_name, path_to_bwa,
                              rg_id, rg_sample, rg_lib, rg_platform):
    '''
    time /exds/users/oabalbin/sw/bwa/bwa-0.5.8a/bwa sampe 
    /exds/projects/alignment_indexes/bwa/hg19_illumina/hg19_illumina.fa 
    ./s_4_1_sequence.txt.sai ./s_4_2_sequence.txt.sai ./s_4_1_sequence.txt 
    ./s_4_2_sequence.txt > s_4_12_sequence.aln.sam
    
    # Using the path for the gatk well formed bam files with read group, 
    and platform information
    time /exds/users/oabalbin/sw/bwa/bwa-0.5.8a/bwa sampe 
    -i MARK -m DEPRISTO -l LIBRARY -p ILLUMINA 
    /exds/projects/alignment_indexes/gatk/hg19/bwa/hg19.fa 
    ./s_3_1_sequence.hg19.sai ./s_3_2_sequence.hg19.sai 
    ./s_3_1_sequence.txt ./s_3_2_sequence.txt > s_3_12_sequence.hg19.aln2.sam
    
    -n INT     Maximum number of alignments to output in the XA tag for reads 
               paired properly. If a read has more than INT hits, the XA tag 
               will not be written. [3]
    -N INT     Maximum number of alignments to output in the XA tag for disconcordant 
               read pairs (excluding singletons). If a read has more than INT hits, 
               the XA tag will not be written. [10] 
    '''
    num_aligns=str(1)
    num_discordant_aligns=str(1)
    #'-n',num_aligns, '-N', num_discordant_aligns,
    #fields=align_reads_pair1.split('_')BWA_read_group_patch.diff
    #name = fields[0]+'_'+fields[1]+'_'+'12'+fields[3]
    #outfile_name=name.replace('.sai','.sam')
    # read group
    read_group="\'@RG\\tID:%s\\tPL:%s\\tLB:%s\\tSM:%s\'"%(rg_id,rg_platform, rg_lib,rg_sample)
    bwa_command = os.path.join(path_to_bwa, 'bwa')
    '''
    args=[bwa_command, 'sampe', '-f', outfile_name, 
          '-i', rg_id, '-m', rg_sample, '-l', rg_lib, '-p', rg_platform,
          ref_genome, align_reads_pair1,align_reads_pair2, raw_reads1, raw_reads2]
    '''
    args=[bwa_command, 'sampe', '-f', outfile_name,
          '-r',read_group, ref_genome, 
          align_reads_pair1,align_reads_pair2, raw_reads1, raw_reads2]
    
    return args
    #args= [a.replace(',',';') for a in args]
    #command = ",".join(args).replace(',',' ').replace(';',',')
    #return command


def main_bwa_alignment(job, config, num_processors, runfunc, deps=None):
    '''
    BWA alignment tasks
    (1) BWA alignment (paired-end or single-read)
    (2) Mark duplicates
    (3) Sort reads in coordinate order
    (4) Index BAM file
    
    job: ExomeLaneConfig object
    config: ExomePipelineConfig object
    num_processors: num of available processors for the job to use
    runfunc: function to use to submit jobs

    returns: tuple (job_id, bam_file) for this job
    '''
    # find reference genome file for species
    genome_tuple = config.genomes[job.species]
    ref_genome = genome_tuple.bwa_ref_genome
    # bwa run parameters
    path_to_bwa = config.bwa_path
    num_mismatches = int(config.bwa_num_mismatches)
    num_cores = int(config.bwa_num_cores)
    # samtools parameters
    path_to_samtools = config.samtools_path
    path_to_picard = config.picard_path
    # email addresses
    email_addresses = config.email_addresses
    # output file name(s)
    
    # TODO: these are 'global' variables and should either be put
    # in a separate config file or moved into the 'base' module    
    node_memory = 45000.0
    node_processors = 12
    single_processor=1
    mem_per_processors = int(float(node_memory) / node_processors)
    extra_mem=8000
    # wt=walltime
    wt_aligning="60:00:00"
    wt_samtools="24:00:00"
    #
    # call the 'bwa aln' command to align each fastq file and create '.sai' files
    # 
    aln_sai_files = []
    aln_job_ids = []
    # TODO: check dependencies before running
    for mate,fastq_file in enumerate(job.fastq_files):
        logging.info('%s: calling bwa aln for mate %d fastq file %s' % (job.name, mate, fastq_file))
        outfile_name = os.path.join(job.align_dir, 'mate%d.sai' % (mate))
        
        if up_to_date(outfile_name, fastq_file):
            logging.info("[SKIPPED] aln file for mate%d.sai is up to date" % (mate))
            aln_sai_files.append(outfile_name)
        else:
            args = reads_alignment(ref_genome, fastq_file, num_mismatches, num_cores, outfile_name, path_to_bwa)
            # TODO: need to specify 'pmem'?
            job_id = runfunc(job.name + '_aln', args, num_cores, cwd=job.output_dir, walltime=wt_samtools, 
                             stdout="bwa_aln_mate%d.log" % (mate),
                             deps=deps, 
                             email_addresses=email_addresses)
            aln_sai_files.append(outfile_name)
            aln_job_ids.append(job_id)
    #
    # call 'bwa samse or sampe' to produce a SAM file
    #
    if all(up_to_date(job.align_sam_file, aln_sai) for aln_sai in aln_sai_files):
        logging.info("[SKIPPED] aln sam file is up to date, %s" % (job.align_sam_file))
    else:
        if len(aln_sai_files) == 1:
            args = single_read_alignment_sam(ref_genome, 
                                             aln_sai_files[0],
                                             job.fastq_files[0],
                                             job.align_sam_file,
                                             path_to_bwa,
                                             rg_id=job.run_id,
                                             rg_sample=job.sample_id,
                                             rg_lib=job.lib_id,
                                             rg_platform=job.seq_platform)
        elif len(aln_sai_files) == 2:
            args = paired_read_alignment_sam(ref_genome, 
                                             aln_sai_files[0], aln_sai_files[1],
                                             job.fastq_files[0], job.fastq_files[1],
                                             job.align_sam_file,
                                             path_to_bwa,
                                             rg_id=job.run_id,
                                             rg_sample=job.sample_id,
                                             rg_lib=job.lib_id,
                                             rg_platform=job.seq_platform)                                         
        else:
            assert False
        # TODO: check dependencies before running
        sam_job_id = runfunc(job.name + "_sam", args, single_processor, cwd=job.output_dir, 
                             walltime=wt_aligning, pmem=extra_mem, 
                             deps=aln_job_ids, 
                             stdout="bwa_sam.log", 
                             email_addresses=email_addresses)
        deps = [sam_job_id]    
    #
    # convert SAM -> BAM
    #
    py_script = os.path.join(_jobs_dir, "sam_to_bam.py")
    args = [sys.executable, py_script, 
            "--sort-order", "coordinate",
            "--samtools-path",path_to_samtools,
            "--picard-path", path_to_picard,
            job.align_sam_file, job.align_bam_file] 
    logging.info("%s: Running SAM -> BAM conversion" % (job.name))
    job_id = runfunc(job.name + "_sb", args, single_processor, cwd=job.output_dir, 
                     walltime=wt_samtools,
                     pmem=extra_mem,
                     stdout="sam2bam.log", 
                     deps=deps, 
                     email_addresses=email_addresses)
    return job_id  
#def bwa_job(lane_file, config_file, num_processors, runfunc, deps=None):



def main_sam2bam_alignment(job, config, num_processors, runfunc, deps=None):
    '''
    BWA alignment tasks
    (1) BWA alignment (paired-end or single-read)
    (2) Mark duplicates
    (3) Sort reads in coordinate order
    (4) Index BAM file
    
    job: ExomeLaneConfig object
    config: ExomePipelineConfig object
    num_processors: num of available processors for the job to use
    runfunc: function to use to submit jobs

    returns: tuple (job_id, bam_file) for this job
    '''
    # find reference genome file for species
    genome_tuple = config.genomes[job.species]
    ref_genome = genome_tuple.bwa_ref_genome
    # bwa run parameters
    path_to_bwa = config.bwa_path
    num_mismatches = int(config.bwa_num_mismatches)
    num_cores = int(config.bwa_num_cores)
    # samtools parameters
    path_to_samtools = config.samtools_path
    path_to_picard = config.picard_path
    # email addresses
    email_addresses = config.email_addresses
    # output file name(s)
    
    # TODO: these are 'global' variables and should either be put
    # in a separate config file or moved into the 'base' module    
    node_memory = 45000.0
    node_processors = 12
    single_processor=1
    mem_per_processors = int(float(node_memory) / node_processors)
    extra_mem=8000
    # wt=walltime
    wt_aligning="60:00:00"
    wt_samtools="24:00:00"
    #
    # call the 'bwa aln' command to align each fastq file and create '.sai' files
    # 
    aln_sai_files = []
    aln_job_ids = []
    # TODO: check dependencies before running
    for mate,fastq_file in enumerate(job.fastq_files):
        logging.info('%s: calling bwa aln for mate %d fastq file %s' % (job.name, mate, fastq_file))
        outfile_name = os.path.join(job.align_dir, 'mate%d.sai' % (mate))
        
        if up_to_date(outfile_name, fastq_file):
            logging.info("[SKIPPED] aln file for mate%d.sai is up to date" % (mate))
            aln_sai_files.append(outfile_name)
        else:
            args = reads_alignment(ref_genome, fastq_file, num_mismatches, num_cores, outfile_name, path_to_bwa)
            # TODO: need to specify 'pmem'?
            job_id = runfunc(job.name + '_aln', args, num_cores, cwd=job.output_dir, walltime=wt_samtools, 
                             stdout="bwa_aln_mate%d.log" % (mate),
                             deps=deps, 
                             email_addresses=email_addresses)
            aln_sai_files.append(outfile_name)
            aln_job_ids.append(job_id)
    #
    # call 'bwa samse or sampe' to produce a SAM file
    #
    if all(up_to_date(job.align_sam_file, aln_sai) for aln_sai in aln_sai_files):
        logging.info("[SKIPPED] aln sam file is up to date, %s" % (job.align_sam_file))
    else:
        if len(aln_sai_files) == 1:
            args = single_read_alignment_sam(ref_genome, 
                                             aln_sai_files[0],
                                             job.fastq_files[0],
                                             job.align_sam_file,
                                             path_to_bwa,
                                             rg_id=job.run_id,
                                             rg_sample=job.sample_id,
                                             rg_lib=job.lib_id,
                                             rg_platform=job.seq_platform)
        elif len(aln_sai_files) == 2:
            args = paired_read_alignment_sam(ref_genome, 
                                             aln_sai_files[0], aln_sai_files[1],
                                             job.fastq_files[0], job.fastq_files[1],
                                             job.align_sam_file,
                                             path_to_bwa,
                                             rg_id=job.run_id,
                                             rg_sample=job.sample_id,
                                             rg_lib=job.lib_id,
                                             rg_platform=job.seq_platform)                                         
        else:
            assert False
        # TODO: check dependencies before running
        sam_job_id = runfunc(job.name + "_sam", args, single_processor, cwd=job.output_dir, 
                             walltime=wt_aligning, pmem=extra_mem, 
                             deps=aln_job_ids, 
                             stdout="bwa_sam.log", 
                             email_addresses=email_addresses)
        deps = [sam_job_id]    
    #
    # convert SAM -> BAM
    #
    py_script = os.path.join(_jobs_dir, "sam_to_bam.py")
    args = [sys.executable, py_script, 
            "--sort-order", "coordinate",
            "--samtools-path",path_to_samtools,
            "--picard-path", path_to_picard,
            job.align_sam_file, job.align_bam_file] 
    logging.info("%s: Running SAM -> BAM conversion" % (job.name))
    job_id = runfunc(job.name + "_sb", args, single_processor, cwd=job.output_dir, 
                     walltime=wt_samtools,
                     pmem=extra_mem,
                     stdout="sam2bam.log", 
                     deps=deps, 
                     email_addresses=email_addresses)
    return job_id


def bwa_job(lane_file, lane_root_dir, config_file, num_processors, runfunc, deps=None):
    config = ExomePipelineConfig()
    config.from_xml(config_file)
    job  = ExomeLaneConfig()
    #job.from_xml(lane_file, config.output_dir)
    job.from_xml(lane_file, lane_root_dir)
    if not os.path.exists(job.output_dir):
        logging.error("%s: output directory %s does not exist" % (job.name, job.output_dir))
        return JOB_ERROR
    # create alignment directory
    if not os.path.exists(job.align_dir):
        logging.info("%s: created alignment directory %s" % (job.name, job.align_dir))
        os.makedirs(job.align_dir)
    # run job
    job_id = main_bwa_alignment(job, config, num_processors, runfunc, deps)
    return job_id
