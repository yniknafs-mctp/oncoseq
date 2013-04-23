'''
Created on Dec 2, 2012

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import xml.etree.cElementTree as etree

from oncoseq.rnaseq.lib.libtable import Library, \
    FRAGMENT_LAYOUT_SINGLE, FRAGMENT_LAYOUT_PAIRED, \
    POLYA_TRANSCRIPTOME, FR_UNSTRANDED
from oncoseq.rnaseq.lib.base import indent_xml

#Code    Definition    Short Letter Code
#01    Primary solid Tumor    TP
#02    Recurrent Solid Tumor    TR
#03    Primary Blood Derived Cancer - Peripheral Blood    TB
#04    Recurrent Blood Derived Cancer - Bone Marrow    TRBM
#05    Additional - New Primary    TAP
#06    Metastatic    TM
#07    Additional Metastatic    TAM
#08    Human Tumor Original Cells    THOC
#09    Primary Blood Derived Cancer - Bone Marrow    TBM
#10    Blood Derived Normal    NB
#11    Solid Tissue Normal    NT
#12    Buccal Cell Normal    NBC
#13    EBV Immortalized Normal    NEBV
#14    Bone Marrow Normal    NBM
#20    Control Analyte    CELLC
#40    Recurrent Blood Derived Cancer - Peripheral Blood    TRB
#50    Cell Lines    CELL
#60    Primary Xenograft Tissue    XP
#61    Cell Line Derived Xenograft Tissue    XCL
SAMPLE_TYPE_MAP = {'01': {'sample_type': 'tissue',
                         'cancer_progression': 'cancer',
                         'recurrence': 'no',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'primary_tumor'},
                  '02': {'sample_type': 'tissue',
                         'cancer_progression': 'cancer',
                         'recurrence': 'yes',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'primary_tumor'},
                  '03': {'sample_type': 'blood',
                         'cancer_progression': 'cancer',
                         'recurrence': 'no',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'blood'},                
                  '04': {'sample_type': 'tissue',
                         'cancer_progression': 'cancer',
                         'recurrence': 'yes',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'bone_marrow'},                
                  '05': {'sample_type': 'tissue',
                         'cancer_progression': 'cancer',
                         'recurrence': 'no',
                         'hx_past_cancer': 'yes',
                         'harvest_site': 'primary_tumor'}, 
                  '06': {'sample_type': 'tissue',
                         'cancer_progression': 'metastatic',
                         'recurrence': 'no',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'metastatic_site'},
                  '07': {'sample_type': 'tissue',
                         'cancer_progression': 'metastatic',
                         'recurrence': 'no',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'metastatic_site'},                
                  '08': {'sample_type': 'cell_line',
                         'cancer_progression': 'cancer',
                         'recurrence': 'no',
                         'hx_past_cancer': 'no'},
                  '09': {'sample_type': 'tissue',
                         'cancer_progression': 'cancer',
                         'recurrence': 'no',
                         'hx_past_cancer': 'no',
                         'harvest_site': 'bone_marrow'},
                  '10': {'sample_type': 'blood',
                         'cancer_progression': 'normal',
                         'harvest_site': 'blood'},
                  '11': {'sample_type': 'tissue',
                         'cancer_progression': 'normal',
                         'harvest_site': 'benign_adjacent'},
                  '12': {'sample_type': 'tissue',
                         'cancer_progression': 'normal',
                         'harvest_site': 'buccal'},
                  '13': {'sample_type': 'cell_line',
                         'cancer_progression': 'normal'},
                  '14': {'sample_type': 'tissue',
                         'cancer_progression': 'normal',
                         'harvest_site': 'bone_marrow'},
                  '20': {'sample_type': 'control',
                         'cancer_progression': 'normal'},
                  '40': {'sample_type': 'blood',
                         'cancer_progression': 'cancer',
                         'recurrence': 'yes',
                         'harvest_site': 'blood'},
                  '50': {'sample_type': 'cell_line'},
                  '60': {'sample_type': 'xenograft'},
                  '61': {'sample_type': 'xenograft_cell_line'}}

#Study Abbreviation    Study Name
#LAML    Acute Myeloid Leukemia
#ACC    Adrenocortical carcinoma
#BLCA    Bladder Urothelial Carcinoma
#LGG    Brain Lower Grade Glioma
#BRCA    Breast invasive carcinoma
#CESC    Cervical squamous cell carcinoma and endocervical adenocarcinoma
#LCLL    Chronic Lymphocytic Leukemia
#LCML    Chronic Myelogenous Leukemia
#COAD    Colon adenocarcinoma
#CNTL    Controls
#ESCA    Esophageal carcinoma 
#GBM    Glioblastoma multiforme
#HNSC    Head and Neck squamous cell carcinoma
#KICH    Kidney Chromophobe
#KIRC    Kidney renal clear cell carcinoma
#KIRP    Kidney renal papillary cell carcinoma
#LIHC    Liver hepatocellular carcinoma
#LUAD    Lung adenocarcinoma
#LUSC    Lung squamous cell carcinoma
#DLBC    Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
#MESO    Mesothelioma
#MISC    Miscellaneous
#OV    Ovarian serous cystadenocarcinoma
#PAAD    Pancreatic adenocarcinoma
#PCPG    Pheochromocytoma and Paraganglioma
#PRAD    Prostate adenocarcinoma
#READ    Rectum adenocarcinoma
#SARC    Sarcoma
#SKCM    Skin Cutaneous Melanoma
#STAD    Stomach adenocarcinoma
#THCA    Thyroid carcinoma
#UCS    Uterine Carcinosarcoma
#UCEC    Uterine Corpus Endometrioid Carcinoma
DISEASE_MAP = {'LAML': {'cancer_type': 'aml'},
               'ACC': {'cancer_type': 'adrenocortical_carcinoma'},
               'BLCA': {'cancer_type': 'bladder_urothelial_carcinoma'},
               'LGG': {'cancer_type': 'brain_glioma',
                       'cancer_subtype': 'lower_grade'},
               'BRCA': {'cancer_type': 'breast_carcinoma',
                        'cancer_subtype': 'invasive'},
               'CESC': {'cancer_type': 'cervical_carcinoma',
                        'cancer_subtype': 'squamous_cell'},
               'LCLL': {'cancer_type': 'cll'},
               'LCML': {'cancer_type': 'cml'},
               'COAD': {'cancer_type': 'colon_adenocarcinoma'},
               'ESCA': {'cancer_type': 'esophageal_carcinoma'},
               'GBM': {'cancer_type': 'glioblastoma_multiforme'},
               'HNSC': {'cancer_type': 'head_neck_squamous_cell_carcinoma'},
               'KICH': {'cancer_type': 'kidney_chromophobe'},
               'KIRC': {'cancer_type': 'kidney_carcinoma',
                        'cancer_subtype': 'clear_cell'},
               'KIRP': {'cancer_type': 'kidney_carcinoma',
                        'cancer_subtype': 'papillary_cell'},
               'LIHC': {'cancer_type': 'liver_hepatocellular_carcinoma'},
               'LUAD': {'cancer_type': 'lung_adenocarcinoma'},
               'LUSC': {'cancer_type': 'lung_squamous_cell_carcinoma'},
               'DLBC': {'cancer_type': 'diffuse_large_bcell_lymphoma'},
               'MESO': {'cancer_type': 'mesothelioma'},
               'OV': {'cancer_type': 'ovarian_serous_cystadenocarcinoma'},
               'PAAD': {'cancer_type': 'pancreatic_adenocarcinoma'},
               'PCPG': {'cancer_type': 'pheochromocytoma_and_paraganglioma'},
               'PRAD': {'cancer_type': 'prostate_adenocarcinoma'},
               'READ': {'cancer_type': 'rectal_adenocarcinoma'},
               'SARC': {'cancer_type': 'sarcoma'},
               'SKCM': {'cancer_type': 'melanoma'},
               'STAD': {'cancer_type': 'gastric_carcinoma'},
               'THCA': {'cancer_type': 'thyroid_carcinoma'},
               'UCS': {'cancer_type': 'uterine_carcinosarcoma'},
               'UCEC': {'cancer_type': 'uterine_corpus_endometrioid_carcinoma'},
               'CNTL': {}}
               
COHORT_MAP = {'LAML': 'blood',
               'ACC': 'adrenal_gland',
               'BLCA': 'bladder',
               'LGG': 'brain',
               'BRCA': 'breast',
               'CESC': 'cervix',
               'LCLL': 'blood',
               'LCML': 'blood',
               'COAD': 'colon',
               'ESCA': 'esophagus',
               'GBM': 'brain',
               'HNSC': 'head_neck',
               'KICH': 'kidney',
               'KIRC': 'kidney',
               'KIRP': 'kidney',
               'LIHC': 'liver',
               'LUAD': 'lung',
               'LUSC': 'lung',
               'DLBC': 'lymphpoid',
               'MESO': 'mesothelium',
               'OV': 'ovary',
               'PAAD': 'pancreas',
               'PCPG': 'ganglia',
               'PRAD': 'prostate',
               'READ': 'colon',
               'SARC': 'connective_tissue',
               'SKCM': 'skin',
               'STAD': 'gastric',
               'THCA': 'thyroid',
               'UCS': 'uterus',
               'UCEC': 'uterus',
               'CNTL': 'control'}               

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", dest="species", default="human")
    parser.add_argument("--library-type", dest="library_type", default=FR_UNSTRANDED)
    parser.add_argument("--library-protocol", dest="library_protocol", default=POLYA_TRANSCRIPTOME)
    parser.add_argument("--param", dest="param_list", action="append", default=None)
    parser.add_argument("--xml", dest="write_xml", action="store_true", default=False)
    parser.add_argument("--ignore", dest="ignore_file", default=None)
    parser.add_argument("cghub_xml_file")
    parser.add_argument("seq_repo")
    parser.add_argument("seq_repo_dir")
    args = parser.parse_args()
    # check params
    if not os.path.exists(args.cghub_xml_file):
        parser.error("cghub_xml_file %s not found" % (args.cghub_xml_file))
    if not os.path.exists(args.seq_repo_dir):
        parser.error("seq_repo_dir %s not found" % (args.seq_repo_dir))
    if args.ignore_file is not None:
        if not os.path.exists(args.ignore_file):
            parser.error("ignore_file %s not found" % (args.ignore_file))
        ignore_libraries = set([x.strip() for x in open(args.ignore_file)])
    else:
        ignore_libraries = set()
    # parse param list
    default_params = {}
    if args.param_list is not None:
        for param in args.param_list:
            k,v = param.split("=")
            default_params[k] = v
    # read libraries
    libraries = {}
    tree = etree.parse(args.cghub_xml_file)  
    root = tree.getroot()
    ignored = 0
    file_not_found = 0
    wrong_filesize = 0
    no_bam = 0
    redundant = 0
    total_results = 0
    no_bam_fileh = open("bam_not_found.txt", "w")
    for elem in root.findall("Result"):
        total_results += 1
        analysis_id = elem.findtext("analysis_id")
        aliquot_id = elem.findtext("aliquot_id")
        if aliquot_id in ignore_libraries:
            logging.warning("Ignoring analysis %s aliquot %s" % (analysis_id, aliquot_id))
            ignored += 1
            continue
        study_id = elem.findtext("study")
        patient_id = elem.findtext("participant_id")
        sample_id = elem.findtext("sample_id")
        legacy_sample_id = elem.findtext("legacy_sample_id")
        disease_abbr = elem.findtext("disease_abbr")
        sample_type = elem.findtext('sample_type')
        analysis_id = analysis_id
        paired_elem = elem.find("experiment_xml/EXPERIMENT_SET/EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED")
        if paired_elem is not None:
            fragment_layout = FRAGMENT_LAYOUT_PAIRED
        else:
            fragment_layout = FRAGMENT_LAYOUT_SINGLE
        bam_files = []
        # find bam file
        files_elem = elem.find("files")
        for file_elem in files_elem.findall("file"):
            filename = file_elem.findtext("filename")
            if os.path.splitext(filename)[-1] != ".bam":
                logging.error("File %s not a BAM file" % (filename))
                continue
            correct_filesize = int(file_elem.findtext("filesize"))
            subpath = os.path.join(analysis_id, filename)
            path = os.path.join(args.seq_repo_dir, subpath)
            if not os.path.exists(path):
                logging.error("Analysis %s not found" % (analysis_id))
                file_not_found += 1
                continue
            filesize = os.path.getsize(path)
            if filesize != correct_filesize:
                logging.error("Analysis %s has incorrect filesize" % (analysis_id))
                wrong_filesize += 1
            else:
                bam_files.append(subpath)
        if len(bam_files) == 0:
            logging.error("Analysis %s has no valid bam files" % (analysis_id))
            no_bam += 1
            print >>no_bam_fileh, analysis_id
            continue
        kwargs = {'study_id': study_id,
                  'cohort_id': COHORT_MAP[disease_abbr],
                  'patient_id': patient_id,
                  'sample_id': sample_id,
                  'library_id': aliquot_id,
                  'description': legacy_sample_id,
                  'species': args.species,
                  'library_protocol': args.library_protocol,
                  'library_type': args.library_type,
                  'fragment_layout': fragment_layout,
                  'seq_repo': args.seq_repo,
                  'read1_files': '',
                  'read2_files': '',
                  'bam_files': ','.join(bam_files)}
        kwargs['params'] = {'tcga_cancer_type': disease_abbr,
                            'tcga_legacy_sample_id': legacy_sample_id,
                            'tcga_sample_type': sample_type}
        kwargs['params'].update(default_params)
        kwargs['params'].update(SAMPLE_TYPE_MAP[sample_type])
        kwargs['params'].update(DISEASE_MAP[disease_abbr])
        library = Library(**kwargs)
        # handle duplicate libraries
        if library.library_id in libraries:
            logging.warning("Redundant result for aliquot %s" % (library.library_id))
            redundant += 1
            # always keep paired end if available
            if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
                libraries[library.library_id] = library
        libraries[library.library_id] = library
    libraries = libraries.values()
    no_bam_fileh.close()
    logging.info("ignored: %d" % (ignored))
    logging.info("file not found: %d" % file_not_found)
    logging.info("wrong filesize: %d" % wrong_filesize)
    logging.info("no bam: %d" % no_bam)
    logging.info("redundant: %d" % (redundant))
    logging.info("total: %d " % (total_results))
    # write
    if args.write_xml:
        root = etree.Element("libraries")
        for library in libraries:
            library.to_xml(root)
        # indent for pretty printing
        indent_xml(root)
        print etree.tostring(root)
    else:
        # build list of all parameters
        params = set()
        for library in libraries:
            params.update(library.params.keys())
        sorted_params = sorted(params)
        # output table
        header_fields = []
        header_fields.extend(Library.fields)
        header_fields.extend(sorted_params)
        print '\t'.join(header_fields)
        for library in libraries:
            fields = []
            library.read1_files = ','.join(library.read1_files)
            library.read2_files = ','.join(library.read2_files)
            library.bam_files = ','.join(library.bam_files)
            for field_name in Library.fields:
                fields.append(getattr(library, field_name))
            for param in sorted_params:
                fields.append(library.params.get(param, "na"))
            print '\t'.join(fields)


if __name__ == '__main__':
    sys.exit(main())
        