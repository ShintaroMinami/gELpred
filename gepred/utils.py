import os
dirpath = os.path.abspath(os.path.dirname(__file__))

FILE_KDINDEX = dirpath+'/data/parameters/Kyte_Doolittle_Hydropathy.param'
DIR_CLREP = dirpath+'/data/clrep/'
DIR_SUBSTR = dirpath+'/data/substr_pdbs/'

MODEL_RW = {
    'RW182': {'pdb': DIR_SUBSTR+'RW182.pdb', 'clrep': DIR_CLREP+'RW182_0.65.clrep', 'kdtid':'KDT1', 'kdcut': 0.65, 'id': 1},
    'RW033': {'pdb': DIR_SUBSTR+'RW033.pdb', 'clrep': DIR_CLREP+'RW033_0.65.clrep', 'kdtid':'KDT2', 'kdcut': 0.65, 'id': 3},
    'RW031': {'pdb': DIR_SUBSTR+'RW031.pdb', 'clrep': DIR_CLREP+'RW031_0.65.clrep', 'kdtid':'KDT3', 'kdcut': 0.65, 'id': 1},
    'RW020': {'pdb': DIR_SUBSTR+'RW020.pdb', 'clrep': DIR_CLREP+'RW020_0.65.clrep', 'kdtid':'KDT4', 'kdcut': 0.65, 'id': 1},
}

MODEL_RWRR = {
    'RW005': {'pdb': DIR_SUBSTR+'RW005.pdb', 'clrep': DIR_CLREP+'RW005_0.60.clrep', 'kdtid':'KDT5', 'kdcut': 0.60, 'id': 1},
    'RR070': {'pdb': DIR_SUBSTR+'RR070.pdb', 'clrep': DIR_CLREP+'RR070_0.60.clrep', 'kdtid':'KDT6', 'kdcut': 0.60, 'id': 3},
    'RR088': {'pdb': DIR_SUBSTR+'RR088.pdb', 'clrep': DIR_CLREP+'RR088_0.60.clrep', 'kdtid':'KDT7', 'kdcut': 0.60, 'id': 1},
    'RR276': {'pdb': DIR_SUBSTR+'RR276.pdb', 'clrep': DIR_CLREP+'RR276_0.60.clrep', 'kdtid':'KDT8', 'kdcut': 0.60, 'id': 1},
    'RW007': {'pdb': DIR_SUBSTR+'RW007.pdb', 'clrep': DIR_CLREP+'RW007_0.60.clrep', 'kdtid':'KDT9', 'kdcut': 0.60, 'id': 2},
}
