'''
Created on Aug 7, 2011
@author: mkiyer
@author: oabalbin
'''
import os
import logging
import xlrd
import xml.etree.cElementTree as etree

# fragment layouts
FRAGMENT_LAYOUT_SINGLE = "single"
FRAGMENT_LAYOUT_PAIRED = "paired"

# quality score formats
SANGER_FORMAT = "sanger"
SOLEXA_FORMAT = "solexa"
ILLUMINA_FORMAT = "illumina"
FASTQ_QUAL_FORMATS = [SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT]

# gender formats
GENDER_MALE = "male"
GENDER_FEMALE = "female"
GENDER_NA = "na"
GENDER_VALUES = [GENDER_MALE, GENDER_FEMALE]

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
    return None

class Patient(object):
    __fields__ = ('id', 'description', 'species', 'ethnicity', 
                  'gender', 'age') 

    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])        
        # relationships
        self.samples = []

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in Patient.__fields__:
            kwargs[f] = elem.findtext(f)
        return Patient(**kwargs)
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "patient")
        for f in Patient.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self,f)
        return parent

    def is_valid(self):
        if self.gender not in GENDER_VALUES:
            logging.error("Patient %s gender %s invalid" % (self.id, self.gender))


class Sample(object):
    __fields__ = ("patient_id", "id", 'description',
                  'sample_type', 'cohort', 'disease',
                  'cancer_progression','protocol')
        
    def __init__(self, **kwargs):
        for attrname in self.__fields__:
            if attrname in kwargs:
                setattr(self, attrname, kwargs[attrname])        
        # relationships
        self.patient = None
        self.libraries = []

    @staticmethod
    def from_xml(elem):
        kwargs = {}
        for f in Sample.__fields__:
            kwargs[f] = elem.findtext(f)
        return Sample(**kwargs)
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "sample")
        for f in Sample.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self,f)
        return parent

    def is_valid(self):
        return True
    


class Library(object):
    __fields__ = ('sample_id', 'id', 'description',
                  'strand_protocol', 'fragment_length','protocol')
        
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
        return Library(**kwargs)
        
    def to_xml(self, root):
        parent = etree.SubElement(root, "library")
        for f in Library.__fields__:
            elem = etree.SubElement(parent, f)
            elem.text = getattr(self, f)
        return parent

    def is_valid(self):
        return True


class Lane(object):    
    __fields__ = ('center_name', 'run', 'lane', 'library_id', 'id',
                  'platform', 'fragment_layout', 'quality_scores',
                  'read1_file', 'read2_file', 'comments', 'qc','exome_kit')

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
        if self.qc == "FAIL":
            logging.error("Lane %s marked QC FAIL" % (self.id))
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
        if not "samples" in sheet_names:
            raise SeqDBError("XLS file missing 'samples' Sheet")
        if not "libraries" in sheet_names:
            raise SeqDBError("XLS file missing 'libraries' Sheet")
        if not "lanes" in sheet_names:
            raise SeqDBError("XLS file missing 'lanes' Sheet")
        # read patients
        patients = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("patients")):
            p = Patient(**field_dict)
            patients[p.id] = p
        # read samples
        samples = {}
        for field_dict in read_wksheet(wkbook.sheet_by_name("samples")):
            #print field_dict
            s = Sample(**field_dict)
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
        # make seqdb object
        seqdb = SeqDB()
        seqdb.patients = patients
        seqdb.samples = samples
        seqdb.libraries = libraries
        seqdb.lanes = lanes
        return seqdb
