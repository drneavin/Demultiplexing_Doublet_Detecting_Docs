#!/usr/env/bin python
import time
import os
import scipy.io
import gzip
import pandas as pd

def import_cellranger_mtx(mtx_directory):
    """
    Reads in a sparse matrix (.mtx) exported by CellRanger in barcodes-by-genes format

    Recognizes output from CellRanger version 2 (files: matrix.mtx) and
    CellRanger v3 (files: matrix.mtx.gz)
    """
    start = time.time()
    if '.mtx' in mtx_directory:
        mtx_file = mtx_directory ### If the file path was all the way to the .mtx file
        mtx_directory = os.path.abspath(os.path.join(mtx_file, os.pardir))
    else:
        mtx_file = os.path.join(mtx_directory, "matrix.mtx")
        
    if not os.path.exists(mtx_file):
        mtx_file = mtx_file + ".gz"
        if not os.path.exists(mtx_file):
            raise Exception("Directory {} does not contain a recognizable matrix file".format(mtx_directory))
    sparse_matrix = scipy.io.mmread(mtx_file).T.tocsc()
    return(sparse_matrix)
    print('sparse matrix imported from mtx file in %s seconds' % str(time.time()-start))


def read_barcodes(barcodes):
    """
    Reads in a sparse matrix (.mtx) exported by CellRanger in barcodes-by-genes format

    Recognizes output from CellRanger version 2 (files: barcodes.tsv) and
    CellRanger v3 (files: barcodes.tsv.gz)
    """
    barcodes_df = pd.read_csv(barcodes, sep = "\t", header=None, names = ["Barcode"])
    return(barcodes_df)