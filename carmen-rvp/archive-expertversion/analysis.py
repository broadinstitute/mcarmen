# General modules
import pandas as pd
import numpy as np
from itertools import groupby
from functools import reduce

def ReadData(x):
    """
    Specifies the rows of the rawdata file x and extracts the raw data for a 192.192 IFC.
    """

    # Set the index within read_csv
    gaps = 2
    top_header = 13
    datasets = 4
    rows = 4608

    # Autocalculating header size here, so easier to fix if format changes
    headers = [top_header + (rows+gaps)*(i+1) for i in range((datasets))]
    
    xprobe = pd.read_csv(x, header = headers[1], nrows = rows, index_col = "Chamber ID")
    xref = pd.read_csv(x, header = headers[0], nrows = rows, index_col = "Chamber ID")
    xprobeb = pd.read_csv(x, header = headers[3], nrows = rows, index_col = "Chamber ID")
    xrefb = pd.read_csv(x, header = headers[2], nrows = rows, index_col = "Chamber ID")
        
    # clean up the data
    for df in [xprobe, xref, xprobeb, xrefb]:
        _drop_lastcolumn(df)
        _format_columns(df)
        
    return xprobe, xref, xprobeb, xrefb


def _drop_lastcolumn(x):
    """
    drops the last column of a dataframe
    """
    x.drop(x.columns[-1], axis=1, inplace=True)
    
def _format_columns(x):
    """
    renames the column x by adding a t
    """
    x.columns = ['t' + str(col).lstrip() for col in x.columns]
    
def LabelInputs(x, y):
    # with x as the signal dataframe
    x.reset_index(inplace=True)
    splitassignment = x['Chamber ID'].str.split("-", n=1, expand=True)
    x["sampleID"] = splitassignment[0]
    x["assayID"] = splitassignment[1]
    x.set_index('Chamber ID', inplace=True)
    # assigns labels from layout sheet y to df x
    s_layout =  pd.read_excel(y, sheet_name='layout_samples', dtype=str)
    a_layout =  pd.read_excel(y, sheet_name='layout_assays', dtype=str)
    # create dictionary with assay or sample numbers and their name
    assays = pd.read_excel(y, sheet_name='assays')
    samples = pd.read_excel(y, sheet_name='samples')
    a_dict = dict(zip(a_layout.values.reshape(-1), assays.values.reshape(-1)))
    s_dict = dict(zip(s_layout.values.reshape(-1), samples.values.reshape(-1)))
    # map
    x['assay'] = x['assayID'].map(a_dict)
    x['sample'] = x['sampleID'].map(s_dict)

def ItemsToList(x, sheet, sort = 'original'):
    y = pd.read_excel(x, sheet_name = sheet)
    if sheet == 'assays':
        new_array = np.stack(y[['C1', 'C2', 'C3']].values, axis=-1)
    if sheet == 'samples':
        new_array = np.stack(y[['C1', 'C2', 'C3', 'C4', 'C5', 'C6', \
                                'C7', 'C8', 'C9', 'C10', 'C11', 'C12', \
                                'C13', 'C14', 'C15', 'C16', 'C17', 'C18', \
                                'C19', 'C20', 'C21', 'C22', 'C23', 'C24']].values, axis=-1)
    itemlist =  np.concatenate(new_array).tolist()
    if sort == 'alphabetical':
        itemlist = np.unique(itemlist)
    if sort == 'original':
        #itemlist = [x[0] for x in groupby(itemlist)]
        itemlist = reduce(lambda l, x: l+[x] if x not in l else l, itemlist, [])
    
    return itemlist

