'''
Created on Aug 7, 2011
@author: mkiyer
@author: oabalbin
'''
import os
import logging
import xlrd
import xml.etree.cElementTree as etree
import collections
import sys

# fragment layouts
FRAGMENT_LAYOUT_SINGLE = "single"
FRAGMENT_LAYOUT_PAIRED = "paired"

# quality score formats
SANGER_FORMAT = "sanger"
SOLEXA_FORMAT = "solexa"
ILLUMINA_FORMAT = "illumina"
FASTQ_QUAL_FORMATS = (SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT)

# source nucleotide types
NUCLEOTIDE_RNA = "rna"
NUCLEOTIDE_DNA = "dna"
NUCLEOTIDE_TYPES = (NUCLEOTIDE_RNA, NUCLEOTIDE_DNA)

# progression types
PROGRESSION_BENIGN = "benign"
PROGRESSION_CANCER = "cancer"
PROGRESSION_METASTATIC = "metastatic"
PROGRESSION_VALUES = (PROGRESSION_BENIGN, PROGRESSION_CANCER, PROGRESSION_METASTATIC)

# strand protocol values
STRAND_PROTOCOL_DUTP = "dutp"
STRAND_PROTOCOL_UNSTRANDED = "unstranded"
STRAND_PROTOCOL_VALUES = (STRAND_PROTOCOL_DUTP, STRAND_PROTOCOL_UNSTRANDED)

# sample group types
SAMPLE_TYPE_EXOME_TUMOR = "exome_tumor"
SAMPLE_TYPE_EXOME_NORMAL = "exome_normal"
SAMPLE_TYPE_RNASEQ = "rnaseq"
SAMPLE_TYPE_CAPTURE_RNASEQ = "capture_rnaseq"
SAMPLE_TYPES = (SAMPLE_TYPE_EXOME_TUMOR, SAMPLE_TYPE_EXOME_NORMAL,
                SAMPLE_TYPE_RNASEQ, SAMPLE_TYPE_CAPTURE_RNASEQ)

# analysis protocol values
PROTOCOL_EXOME_DNA = "exome"
PROTOCOL_POLYA_RNA = "rnaseq"
PROTOCOL_EXOME_RNA = "capture_rnaseq"
VALID_PROTOCOLS = (PROTOCOL_EXOME_DNA, PROTOCOL_POLYA_RNA, PROTOCOL_EXOME_RNA)

# Exome kits
VALID_EXOME_KITS={'agilent_v4':"capture_agilent", 'agilent_v3':"capture_agilent_v3",'roche_v2':"capture_roche", 'truseq':"capture_truseq"}


class SeqDBError(Exception):
    pass

def _find_sequence_file(filename):
    if os.path.isfile(filename):
        return filename
    newfile = filename + ".gz"
    if os.path.isfile(newfile):
        return newfile
    newfile = os.path.splitext(filename)[0]
    if os.path.isfile(newfile):
        return newfile
    else:
        # Dangerous
        newfile = filename.split("/")[-1]
        return newfile
    return None

class Patient(object):
    __fields__ = ('id', 'description', 'species', 'study')

    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])        
        # custom parameters
        self.params = kwargs["params"]
        # relationships
        self.samples = []
        self.sample_groups = {}

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in Patient.__fields__:
            kwargs[f] = elem.findtext(f)
        params = {}
        for param_elem in elem.findall("param"):
            params[param_elem.get("name")] = param_elem.text
        kwargs["params"] = params
        patient = Patient(**kwargs)
        sample_dict = {}
        # read samples
        for sample_elem in elem.findall("sample"):
            sample = Sample.from_xml(sample_elem)
            sample_dict[sample.id] = sample
            # link to patient
            sample.patient = patient
            patient.samples.append(sample)
        # read sample groups
        for grp_elem in elem.findall("sample_group"):
            grp = SampleGroup.from_xml(grp_elem)
            for sample_type in SAMPLE_TYPES:
                v = getattr(grp, sample_type)
                if not v:
                    setattr(grp, sample_type, None)
                    grp.samples[sample_type] = None
                elif v not in sample_dict:
                    logging.error("SampleGroup %s cannot link to %s sample %s" % (grp.id, sample_type, v))
                    setattr(grp, sample_type, None)
                    grp.samples[sample_type] = None
                else:
                    grp.samples[sample_type] = sample_dict[v]
            grp.patient = patient
            patient.sample_groups[grp.id] = grp
        return patient

    def to_xml(self, root):
        parent = etree.SubElement(root, "patient")
        for f in Patient.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self,f)
        for k,v in self.params.iteritems():
            elem = etree.SubElement(parent, "param", name=k)
            elem.text = v
        for grp in self.sample_groups.itervalues():
            grp.to_xml(parent)
        # recurse
        for sample in self.samples:
            sample.to_xml(parent)
        return parent

    def is_valid(self):
        is_valid = True
        for sample in self.samples:
            is_valid = is_valid and sample.is_valid()
        return is_valid


class Sample(object):
    __fields__ = ("patient_id", "id", 'description',
                  'sample_type', 'nucleotide_type',
                  'cohort', 'disease',
                  'cancer_progression')

    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])
        
        # custom parameters
        self.params = kwargs["params"]
        # relationships
        self.patient = None
        self.libraries = []
        

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in Sample.__fields__:
            kwargs[f] = elem.findtext(f)
        params = {}
        for param_elem in elem.findall("param"):
            params[param_elem.get("name")] = param_elem.text
        kwargs["params"] = params
        sample = Sample(**kwargs)
        # read libraries
        for lib_elem in elem.findall("library"):
            lib = Library.from_xml(lib_elem)
            lib.sample = sample
            sample.libraries.append(lib)
        return sample
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "sample")
        for f in Sample.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self,f)
        for k,v in self.params.iteritems():
            elem = etree.SubElement(parent, "param", name=k)
            elem.text = v
        # recurse
        for library in self.libraries:
            library.to_xml(parent)
        return parent

    def is_valid(self):
        is_valid = True
        if self.nucleotide_type not in NUCLEOTIDE_TYPES:
            logging.error("Invalid nucleotide type %s" % (self.nucleotide_type))
            is_valid = False
        if self.cancer_progression not in PROGRESSION_VALUES:
            logging.error("Invalid progression value %s" % (self.cancer_progression))
            is_valid = False
        if self.patient is None:
            logging.error("Sample %s unknown patient" % (self.id))
            is_valid = False            
        for library in self.libraries:
            is_valid = is_valid and library.is_valid()
        return is_valid

class SampleGroup(object):    
    __fields__ = ('patient_id', 'id', 'exome_tumor', 'exome_normal',
                  'rnaseq', 'capture_rnaseq')

    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])
        # relationships
        self.patient = None
        self.samples = {}

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in SampleGroup.__fields__:
            kwargs[f] = elem.findtext(f)
        return SampleGroup(**kwargs)
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "sample_group")
        for f in SampleGroup.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self, f)
        return parent

    def is_valid(self):
        valid = True
        return valid

class Library(object):
    __fields__ = ('sample_id', 'id', 'description',
                  'strand_protocol', 'fragment_length', 
                  'capture_kit', 'protocol')
        
    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])        
        self.fragment_length = str(float(self.fragment_length))
        # relationships
        self.sample = None
        self.lanes = []

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in Library.__fields__:
            kwargs[f] = elem.findtext(f)
        lib = Library(**kwargs)
        # read lanes
        for lane_elem in elem.findall("lane"):
            lane = Lane.from_xml(lane_elem)
            lane.library = lib
            #print lane.library
            lib.lanes.append(lane)
        return lib
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "library")
        for f in Library.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self, f)
        # recurse
        for lane in self.lanes:
            # TODO: check to make sure this is the best way to exclude 
            # QC fail libraries
            if lane.qc == "FAIL":
                continue
            lane.to_xml(parent)
        return parent

    def is_valid(self):
        is_valid = True
        if self.protocol not in VALID_PROTOCOLS:
            logging.error("Invalid protocol value %s" % (self.protocol))
            is_valid = False
        # TODO: add capture kit check here
        # check that RNA samples have valid strand protocol
        if ((self.protocol != PROTOCOL_EXOME_DNA) and
            self.strand_protocol not in STRAND_PROTOCOL_VALUES):
            logging.error("Invalid strand protocol %s" % (self.strand_protocol))
            is_valid = False
        if self.sample is None:
            logging.error("Library %s unknown sample" % (self.id))
            is_valid = False            
        for lane in self.lanes:
            is_valid = is_valid and lane.is_valid()
        return is_valid


class Lane(object):    
    __fields__ = ('center_name', 'run', 'lane', 'library_id', 'id',
                  'platform', 'fragment_layout', 'quality_scores',
                  'read1_file', 'read2_file', 'comments', 'qc')

    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])
        self.lane = self.lane
        self.read1_file = _find_sequence_file(self.read1_file)
        self.read2_file = _find_sequence_file(self.read2_file)
        # relationships
        self.library = None

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in Lane.__fields__:
            kwargs[f] = elem.findtext(f)
        return Lane(**kwargs)
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "lane")
        for f in Lane.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self, f)
        return parent

    def is_valid(self):
        valid = True
        if (self.read1_file is None) or (not os.path.exists(self.read1_file)):
            logging.error("Lane %s read 1 file %s not found" % (self.id, self.read1_file))
            valid = False
        if (self.fragment_layout == FRAGMENT_LAYOUT_PAIRED):
            if ((self.read2_file is None) or (not os.path.exists(self.read2_file))):
                logging.error("Lane %s read 2 file %s not found" % (self.id, self.read2_file))
                valid = False
            elif self.read1_file == self.read2_file:
                logging.error("Lane %s read 1 file %s is same as read 2 file" % (self.id, self.read1_file))
                valid = False                
        if self.quality_scores not in FASTQ_QUAL_FORMATS:
            logging.error("Lane %s unknown quality scores format %s" % (self.id, self.quality_scores))
            valid = False
        if self.library is None:
            logging.error("Lane %s unknown library" % (self.id))
            valid = False
        return valid

def read_wksheet(wksheet):
    field_names = wksheet.row_values(0)
    field_descs = wksheet.row_values(1)
    for rownum in xrange(2, wksheet.nrows):
        fields = wksheet.row_values(rownum)
        #print fields
        fields = [' '.join(str(field).split('\n')) for field in fields]
        # build dictionary of field names to field values
        field_name_value_dict = dict((field_names[i], fields[i]) for i in xrange(len(fields)))
        yield field_name_value_dict

class SeqDB(object):
    @staticmethod
    def from_xls(filename):
        """
        parses an XLS file and constructs 'Patient', 'Sample', 'Library', 
        and 'Lane' objects stored as dictionaries keyed by unique id 
        """
        if not os.path.isfile(filename):
            raise OSError("File %s not found or not a regular file" % (filename))
        wkbook = xlrd.open_workbook(filename)
        # check that required sheet names exist
        sheet_names = wkbook.sheet_names()
        if not "patients" in sheet_names:
            raise SeqDBError("XLS file missing 'patients' Sheet")
        if not "patient_parameters" in sheet_names:
            raise SeqDBError("XLS file missing 'patient_parameters' Sheet")
        if not "samples" in sheet_names:
            raise SeqDBError("XLS file missing 'samples' Sheet")
        if not "sample_parameters" in sheet_names:
            raise SeqDBError("XLS file missing 'sample_parameters' Sheet")
        if not "libraries" in sheet_names:
            raise SeqDBError("XLS file missing 'libraries' Sheet")
        if not "lanes" in sheet_names:
            raise SeqDBError("XLS file missing 'lanes' Sheet")
        # read patient parameters
        patient_params = collections.defaultdict(lambda: {})
        for field_dict in read_wksheet(wkbook.sheet_by_name("patient_parameters")):
            patient_id = field_dict["id"]
            k = field_dict["parameter_name"]
            v = field_dict["parameter_value"]
            patient_params[patient_id][k] = v
        # read patients
        patients = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("patients")):
            # add params to field dict
            field_dict["params"] = patient_params[field_dict["id"]]
            # build patient object
            p = Patient(**field_dict)
            # ensure unique ids
            if p.id in patients:
                #print p.id
                raise SeqDBError("Found duplicate patient id %s" % (p.id))
            patients[p.id] = p
        # read sample parameters
        sample_params = collections.defaultdict(lambda: {})
        for field_dict in read_wksheet(wkbook.sheet_by_name("sample_parameters")):
            sample_id = field_dict["id"]
            k = field_dict["parameter_name"]
            v = field_dict["parameter_value"]
            sample_params[sample_id][k] = v
        # read samples
        samples = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("samples")):
            # add params to field dict
            field_dict["params"] = sample_params[field_dict["id"]]
            # build sample object
            s = Sample(**field_dict)
            # ensure unique ids
            if s.id in samples:
                raise SeqDBError("Found duplicate sample id %s" % (s.id))
            # link samples to patients
            p = patients[s.patient_id]
            p.samples.append(s)
            s.patient = p
            samples[s.id] = s
        # read libraries
        libraries = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("libraries")):
            lib = Library(**field_dict)
            s = samples[lib.sample_id]
            s.libraries.append(lib)
            lib.sample = s
            libraries[lib.id] = lib
        # read lanes        
        lanes = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("lanes")):
            lane = Lane(**field_dict)
            lib = libraries[lane.library_id]
            lib.lanes.append(lane)
            lane.library = lib
            if not lane.is_valid():
                logging.error("Lane %s skipped" % (lane.id))
            lanes[lane.id] = lane
        # read sample groups
        sample_groups = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("groups")):
            grp = SampleGroup(**field_dict)
            p = patients[grp.patient_id]
            p.sample_groups[grp.id] = grp
            grp.patient = p
            for sample_type in SAMPLE_TYPES:
                v = getattr(grp, sample_type)
                if not v:
                    setattr(grp, sample_type, None)
                    grp.samples[sample_type] = None
                elif v not in samples:
                    logging.error("Cannot link to sample %s in SampleGroup %s" % (v, grp.id))
                    setattr(grp, sample_type, None)
                    grp.samples[sample_type] = None
                else:
                    grp.samples[sample_type] = samples[v]
            sample_groups[grp.id] = grp
        # make seqdb object
        seqdb = SeqDB()
        seqdb.patients = patients
        seqdb.sample_groups = sample_groups
        seqdb.samples = samples        
        seqdb.libraries = libraries
        seqdb.lanes = lanes
        return seqdb
