'''
Created on Feb 9, 2012

@author: mkiyer
'''
import numpy as np
import collections

NUM_METACOLS = 5

class ExpressionMatrix(object):    

    def __init__(self):
        self.tags = None
        self.metacol_names = None
        self.metacols = None
        self.colnames = None
        self.metarows = None
        self.data = None
        self.sizefactors = None
        self.scalefunc = None

    @staticmethod        
    def fromtabular(fh, sep="\t"):
        mat = ExpressionMatrix()
        # read header
        header = fh.next().strip()
        header_fields = header.split(sep)
        # set column names
        mat.metacol_names = header_fields[:NUM_METACOLS]
        mat.colnames = header_fields[NUM_METACOLS:]
        # read metadata rows
        mat.metarows = {}
        for line in fh:
            fields = line.rstrip().split(sep)
            if fields[0]:
                break
            mat.metarows[fields[NUM_METACOLS-1]] = fields[NUM_METACOLS:]
        # read first line
        mat.metacols = []
        mat.data = []
        mat.metacols.append(fields[:NUM_METACOLS])
        mat.data.append(map(int,fields[NUM_METACOLS:]))
        for line in fh:
            fields = line.rstrip().split(sep)
            mat.metacols.append(fields[:NUM_METACOLS])
            mat.data.append(map(int,fields[NUM_METACOLS:]))
        # convert data from list to numpy array
        mat.data = np.array(mat.data, dtype=np.int)
        return mat

    def merge_cols(self, attrname, merge_func=np.median):
        # group fields together
        attrvals = self.metarows[attrname]
        col_group_dict = collections.OrderedDict()
        for i,val in enumerate(attrvals):
            if val not in col_group_dict:                
                col_group_dict[val] = []
            col_group_dict[val].append(i)
        
            
        
            
        
        
#        # read tags
#        tagdict = {}
#        for line in fh:
#            if not line.startswith("@"):
#                break
#            tag, val = parse_tag(line)
#            tagdict[tag] = val
#        mat.tags = tagdict
#        # next line must have column names
#        colfields = line.strip().split(sep)
#        # next line must have first metadata row
#        line = fh.next().rstrip()
#        fields = line.split(sep)
#        # find first non-empty field
#        for i,f in enumerate(fields):
#            if f: break
#        num_metacols = i + 1        
#        # now get column names
#        colnames = collections.OrderedDict()
#        for i,col in enumerate(colfields[num_metacols:]):
#            colnames[col] = i
#        mat.colnames = colnames
#        num_datacols = len(colfields) - num_metacols
#        # read metarows into standard python dict
#        metarows = {}
#        while True:                
#            fields = line.rstrip().split(sep)
#            if fields[0]:
#                break
#            name = fields[num_metacols-1]            
#            if len(fields[num_metacols:]) == 0:
#                maxlen = 64
#            else:                
#                maxlen = max(len(f) for f in fields[num_metacols:])             
#            vals = np.zeros(num_datacols, dtype="a%s" % (maxlen))
#            vals[:len(fields[num_metacols:])] = fields[num_metacols:]
#            metarows[name] = vals
#            line = fh.next()
#        mat.metarows = metarows
#        # read data
#        lines = [line]
#        lines.extend(line for line in fh)
#        num_datarows = len(lines)
#        num_datacols = len(lines[0].strip().split(sep)) - num_metacols
#        # allocate record array for meta cols
#        metacol_descr = mat.tags["MCOL_CLASSES"].split(",")
#        mat.metacol_classes = [np.cast[col] for col in metacol_descr]
#        metacol_dtype = zip(colfields[:num_metacols], metacol_descr)
#        mat.metacols = np.zeros(num_datarows, dtype=metacol_dtype)
#        # allocate array for data
#        mat.data = np.zeros((num_datarows, num_datacols), dtype=np.float)
#        # convert to matrix format
#        for rownum,line in enumerate(lines):
#            fields = line.strip().split(sep)
#            rowmetacols = tuple(mat.metacol_classes[i](f) for i,f in enumerate(fields[:num_metacols]))
#            rowdata = map(float, fields[num_metacols:])
#            mat.metacols[rownum] = rowmetacols
#            mat.data[rownum,:] = rowdata
#        return mat

# tracking_id     locus   nearest_ref_id  class_code      transcript_length
