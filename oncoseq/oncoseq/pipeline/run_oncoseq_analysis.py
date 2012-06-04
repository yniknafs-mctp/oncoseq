'''
Created on May 17, 2012

@author: mkiyer
'''
import sys
import os
import argparse
import logging
import xml.etree.cElementTree as etree

from oncoseq.lib import config
from oncoseq.lib.config import AnalysisConfig, PipelineConfig
from oncoseq.lib.cluster import submit_job_pbs, submit_job_nopbs
from oncoseq.lib.base import indent_xml

from oncoseq.rnaseq.pipeline import run_rnaseq_analysis
from oncoseq.exome.pipeline import run_exome_analysis
import oncoseq.pipeline
_pipeline_dir = oncoseq.pipeline.__path__[0]

def run_sample_group(grp, genome, server, pipeline, 
                     num_processors, submit_job_func, keep_tmp):
    logging.info("Analyzing sample group: %s" % (grp.id)) 
    # create output dir
    if not os.path.exists(grp.output_dir):
        logging.info("Creating sample group directory: %s" % (grp.output_dir))
        os.makedirs(grp.output_dir)
    # output sample XML
    root = etree.Element("analysis")
    grp.to_xml(root)
    indent_xml(root)
    f = open(grp.xml_file, "w")
    print >>f, etree.tostring(root)
    f.close()   
    # process rnaseq/exome samples
    deps = []
    deps.extend(run_rnaseq_analysis.run_sample_group(grp, genome, server, 
                                                     pipeline, 
                                                     num_processors, 
                                                     submit_job_func, 
                                                     keep_tmp))
    deps.extend(run_exome_analysis.run_sample_group(grp, genome, server, 
                                                    pipeline, num_processors, 
                                                    submit_job_func, 
                                                    keep_tmp))
    return deps


def run_patient(patient, server, pipeline, num_processors,
                submit_job_func, keep_tmp):
    logging.info("Analyzing patient: %s" % (patient.id)) 
    # check genome
    if patient.species not in pipeline.species:
        logging.error("Patient %s genome %s not found" % 
                      (patient.id, patient.species))
        return config.JOB_ERROR
    # get genome corresponding to patient
    genome = pipeline.species[patient.species]
    # create output dir
    if not os.path.exists(patient.output_dir):
        logging.info("Creating patient directory: %s" % (patient.output_dir))
        os.makedirs(patient.output_dir)
    # output patient XML
    root = etree.Element("analysis")
    patient.to_xml(root)
    indent_xml(root)
    f = open(patient.xml_file, "w")
    print >>f, etree.tostring(root)
    f.close()
    #
    # process patient sample groups
    #
    grp_deps = []
    for grp in patient.sample_groups.itervalues():
        logging.info("Running patient %s sample group %s" % (patient.id, grp.id))
        grp_deps.extend(run_sample_group(grp, genome, server, pipeline, 
                                         num_processors, submit_job_func, keep_tmp))
    #
    # write file indicating patient job is complete
    #
    msg = "Notifying user that patient job is complete"
    deps = []
    if os.path.exists(patient.job_complete_file) and (len(grp_deps) == 0):
        logging.info("[SKIPPED]: %s" % msg)
    else:
        logging.info(msg)
        args = [sys.executable, os.path.join(_pipeline_dir, "notify_complete.py"),
                patient.job_complete_file]
        job_id = submit_job_func("done_%s" % (patient.id), args,
                                 num_processors=1,
                                 node_processors=server.node_processors,
                                 node_memory=server.node_mem,
                                 pbs_script_lines=server.pbs_script_lines,
                                 working_dir=patient.output_dir,
                                 walltime="1:00:00",
                                 email="ae",
                                 deps=grp_deps)
        deps = [job_id]
    return deps

def run_analysis(analysis_file, config_file, server_name,
                 num_processors, keep_tmp):
    """
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
    if not pipeline.is_valid(server_name):
        logging.error("Pipeline config not valid")
        return config.JOB_ERROR
    # setup server
    server = pipeline.servers[server_name]
    if server.pbs:
        submit_job_func = submit_job_pbs
    else:
        submit_job_func = submit_job_nopbs
    #
    # attach analysis to pipeline output directory
    #
    analysis.attach_to_results(server.output_dir)    
    #
    # process each patient in analysis
    #
    for patient in analysis.patients:
        run_patient(patient, server, pipeline, 
                    num_processors, submit_job_func, keep_tmp)
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--keep-tmp", action="store_true", dest="keep_tmp", default=False)
    parser.add_argument("analysis_file")
    parser.add_argument("config_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    return run_analysis(args.analysis_file, args.config_file, 
                        args.server_name, 
                        num_processors=args.num_processors,
                        keep_tmp=args.keep_tmp)
    
if __name__ == '__main__': 
    sys.exit(main())