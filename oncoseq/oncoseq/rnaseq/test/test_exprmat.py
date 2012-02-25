'''
Created on Feb 9, 2012

@author: mkiyer
'''
import sys
import numpy as np
import collections

from oncoseq.lib.exprmat import ExpressionMatrix

def sizefunc(m,i):
    colind = m.metacol_names.index("transcript_length")
    return int(m.metacols[i][colind]) > 200

def totalfragfunc(m,i):
    return int(m.metarows["total_frags"][i]) >= 5000000

def tumorfunc(m,i):
    return m.metarows["progression"][i] == "cancer"

def myfunc(m,i):
    return ((m.metarows["sample_type"][i] != "cell_line") and 
            (m.metarows["progression"][i] != "benign") and
            (not m.metarows["sample"][i].startswith("mo_")))

def reorder_by_metarows(m, rownames, dtypes):
    groups = collections.defaultdict(lambda: [])
    for i in xrange(len(m.colnames)):
        key = tuple(dtypes[x](m.metarows[rownames[x]][i]) for x in xrange(len(rownames)))
        groups[key].append(m.colnames[i])
    reordered_colnames = []
    for k in sorted(groups):
        reordered_colnames.extend(sorted(groups[k]))
    colmap = []
    for colname in reordered_colnames:
        colmap.append(m.colnames.index(colname))
    return colmap

def reorder_by_datarow(m, tracking_id, reverse=False):
    row_index = None
    for i,metacolrow in enumerate(m.metacols):
        if metacolrow[0] == tracking_id:
            row_index = i
            break
    colmap = np.argsort(m.data[row_index,:])
    if reverse:
        colmap = colmap[::-1]
    return colmap

if __name__ == '__main__':
    filename = sys.argv[1]
    mat = ExpressionMatrix.fromtabular(open(filename))
    # remove libraries with few fragments
    newmat = mat.filter_cols(totalfragfunc)
    # get just tumor samples
    newmat = newmat.filter_cols(tumorfunc)
    # get met and tumor tissue samples (not mi-oncoseq)
    #newmat = newmat.filter_cols(myfunc)    
    # remove small genes
    countmat = newmat.filter_rows(sizefunc)
    rpkmmat = newmat.filter_rows(sizefunc)
    # group same samples together    
    countmat = countmat.merge_cols("sample")
    # convert to rpkm and group samples
    rpkmmat = rpkmmat.to_rpkm([1e6*int(x) for x in rpkmmat.metarows["total_frags"]])
    rpkmmat = rpkmmat.merge_cols("sample")
    # get column order
    #col_map = reorder_by_metarows(countmat, ("sample_type","progression"), (str, str))    
    #col_map = reorder_by_datarow(rpkmmat, "G48505")    
    #col_map = reorder_by_datarow(rpkmmat, "G85450")    
    countmat = countmat.reorder_cols(col_map)
    # print matrix
    for line in countmat.to_tabular():
        print line

        
            
