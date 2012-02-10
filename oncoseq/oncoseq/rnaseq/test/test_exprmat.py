'''
Created on Feb 9, 2012

@author: mkiyer
'''
import sys

from oncoseq.lib.exprmat import ExpressionMatrix

if __name__ == '__main__':
    filename = sys.argv[1]
    mat = ExpressionMatrix.fromtabular(open(filename))
    print mat.metacol_names
    print mat.colnames
    print mat.metarows.keys()
    
    import numpy as np
    print np.sum(mat.data, axis=0)