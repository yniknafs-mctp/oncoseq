'''
Created on Feb 9, 2012

@author: mkiyer
'''
import numpy as np
import collections
import copy

NUM_METACOLS = 5

ROW_TRACKING_ID = 0
ROW_LOCUS = 1
ROW_NEAREST_REF_ID = 2
ROW_CLASS_CODE = 3
ROW_TRANSCRIPT_LENGTH = 4

DTYPE_RPKM = "rpkm"
DTYPE_COUNT = "count"
DTYPE_FUNCS = {DTYPE_RPKM: float, DTYPE_COUNT: int}

class ExpressionMatrix(object):    

    def __init__(self):
        self.metacol_names = None
        self.metacols = None
        self.colnames = None
        self.metarows = None
        self.data = None

    @staticmethod        
    def fromtabular(fh, dtype=DTYPE_COUNT, sep="\t"):
        dtype_func = DTYPE_FUNCS[dtype]
        # read header
        header = fh.next().strip()
        header_fields = header.split(sep)
        # set column names
        mat = ExpressionMatrix()
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
        replace_missing = lambda fields: [("-" if not x else x) for x in fields]
        mat.metacols.append(replace_missing(fields[:NUM_METACOLS]))
        mat.data.append(map(dtype_func,fields[NUM_METACOLS:]))
        for line in fh:
            fields = line.rstrip().split(sep)
            mat.metacols.append(replace_missing(fields[:NUM_METACOLS]))
            mat.data.append(map(dtype_func,fields[NUM_METACOLS:]))
        # convert data from list to numpy array
        mat.data = np.array(mat.data, dtype=dtype_func)
        return mat
    
    def to_tabular(self, sep="\t"):
        # header
        yield sep.join(self.metacol_names + self.colnames)
        # metadata rows
        blank_fields = [""] * (NUM_METACOLS-1)
        for attrname,attrvals in self.metarows.iteritems():
            yield sep.join(blank_fields + [attrname] + attrvals)
        # data rows
        for i in xrange(self.data.shape[0]):
            yield sep.join(self.metacols[i] + map(str, self.data[i,:]))

    def filter_cols(self, filter_func):
        # make new matrix 
        mat = ExpressionMatrix()
        mat.metacol_names = list(self.metacol_names)
        mat.metacols = copy.deepcopy(self.metacols)
        col_map = []
        for i in xrange(len(self.colnames)):
            if filter_func(self, i):
                col_map.append(i)
        mat.colnames = [self.colnames[x] for x in col_map]
        mat.metarows = {}
        for tag, vals in self.metarows.iteritems():
            mat.metarows[tag] = [vals[x] for x in col_map]
        # copy data        
        nrows = self.data.shape[0]
        ncols = len(mat.colnames)
        mat.data = np.zeros((nrows,ncols), dtype=self.data.dtype)
        for i in xrange(len(mat.colnames)):
            mat.data[:,i] = self.data[:,col_map[i]]
        return mat

    def filter_rows(self, filter_func):
        # make new matrix 
        mat = ExpressionMatrix()
        mat.metacol_names = list(self.metacol_names)
        mat.colnames = list(self.colnames)
        mat.metarows = copy.copy(self.metarows)
        # find which rows to keep
        keep_rows = []
        for i in xrange(self.data.shape[0]):
            if filter_func(self, i):
                keep_rows.append(i)
        # get filtered matrix
        mat.metacols = []
        mat.data = np.zeros((len(keep_rows),self.data.shape[1]), 
                            dtype=self.data.dtype)
        for i,rowindex in enumerate(keep_rows):
            mat.metacols.append(list(self.metacols[rowindex]))
            mat.data[i,:] = self.data[rowindex,:]
        return mat

    def reorder_cols(self, col_map):
        # make new matrix 
        mat = ExpressionMatrix()
        mat.metacol_names = list(self.metacol_names)
        mat.metacols = copy.deepcopy(self.metacols)
        mat.colnames = [self.colnames[x] for x in col_map]
        mat.metarows = {}
        for tag, vals in self.metarows.iteritems():
            mat.metarows[tag] = [vals[x] for x in col_map]
        # copy data
        nrows = self.data.shape[0]
        ncols = len(col_map)  
        mat.data = np.zeros((nrows,ncols), dtype=self.data.dtype)
        for i in xrange(len(col_map)):
            mat.data[:,i] = self.data[:,col_map[i]]
        return mat

    def merge_cols(self, attrname, merge_func=np.median):
        # group fields together
        attrvals = self.metarows[attrname]
        col_group_dict = collections.OrderedDict()
        for i,val in enumerate(attrvals):
            if val not in col_group_dict:                
                col_group_dict[val] = []
            col_group_dict[val].append(i)
        # make new matrix 
        mat = ExpressionMatrix()
        mat.metacol_names = list(self.metacol_names)
        mat.metacols = copy.deepcopy(self.metacols)
        mat.colnames = list(col_group_dict.keys())
        # merge metadata rows and data values
        mat.metarows = {}
        # add a metadata row for the current column names
        newvals = []
        for colname in mat.colnames:
            indexes = col_group_dict[colname]
            newval = ",".join((self.colnames[i] for i in indexes))
            newvals.append(newval)
        mat.metarows["old_colnames"] = newvals
        for tag, vals in self.metarows.iteritems():
            # skip the attribute that we are merging on
            if tag == attrname:
                continue
            newvals = []
            for colname in mat.colnames:
                indexes = col_group_dict[colname]
                newval = ",".join(sorted(set(vals[i] for i in indexes)))
                #newval = ",".join((vals[i] for i in indexes))
                newvals.append(newval)
            mat.metarows[tag] = newvals
        # merge data        
        nrows = self.data.shape[0]
        ncols = len(mat.colnames)
        mat.data = np.zeros((nrows,ncols),dtype=self.data.dtype)
        for i,colname in enumerate(mat.colnames):
            indexes = np.array(col_group_dict[colname])
            mat.data[:,i] = merge_func(self.data[:,indexes], axis=1)
        return mat

    def to_rpkm(self, sizefactors, exp=3):
        # make a new matrix
        mat = ExpressionMatrix()
        mat.metacol_names = list(self.metacol_names)
        mat.metacols = copy.deepcopy(self.metacols)
        mat.metarows = copy.deepcopy(self.metarows)
        mat.colnames = list(self.colnames)
        mat.data = np.zeros(self.data.shape, dtype=DTYPE_FUNCS[DTYPE_RPKM])
        for i in xrange(len(mat.colnames)):
            mat.data[:,i] = self.data[:,i]
        # prepare a matrix that will normalize counts by library
        # size and feature size (e.g. RPKM)
        sizefactors = np.array(sizefactors, dtype=np.float)
        sizefactors = np.reshape(sizefactors, (1,self.data.shape[1]))
        transcript_lengths = np.array([int(row[ROW_TRANSCRIPT_LENGTH]) for row in self.metacols], dtype=np.float)
        transcript_lengths = np.reshape(transcript_lengths, (self.data.shape[0],1))
        norm_matrix = 1.0**(exp) / (transcript_lengths * sizefactors)
        mat.data *= norm_matrix
        mat.dtype = DTYPE_RPKM
        return mat
