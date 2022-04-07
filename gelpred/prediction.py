import subprocess
import datetime
from tqdm import tqdm
from .utils import MODEL_RW, MODEL_RWRR
from .filtstruct import get_structure_hits
from .filtkdind import filter_by_kdindex
import os

MAINDIR = os.path.abspath(os.path.dirname(__file__))
MICAN = MAINDIR+'/bin/mican_2015.03.09'
MICAN_OPTION = {'RW': '-w -g 100 -n 5', 'RR': '-r -g 100 -n 5'}
MODELS = {'RW': MODEL_RW, 'RWRR': MODEL_RWRR}


def calc_mican(pdbfiles, models=None, tempdir='./tmp/', pbar=None, no_mican_warning=True, return_tmpfiles=False):
    stderr = '2> /dev/null' if no_mican_warning else ''
    modeltypes = ['RW', 'RWRR'] if models == None else models
    mican_result = []
    tmp_files = []
    for modeltype in modeltypes:
        model = MODELS[modeltype]
        substrlist = []
        for key in model.keys():
            alnlist = []
            for pdb in pdbfiles:
                if pbar: pbar.update(1)
                substr = model[key]
                substr_pdb = substr['pdb']
                protname = os.path.splitext(os.path.basename(pdb))[0]
                strdatetime = datetime.datetime.now().strftime('@%Y:%m:%d:%H:%M:%S')
                alnfile = tempdir + '@'.join([key, protname]) + strdatetime +'.aln'
                subprocess.check_output('{} {} {} {} -a {} {}'.format(MICAN, substr_pdb, pdb, MICAN_OPTION[key[0:2]], alnfile, stderr), text=True, shell=True)
                alnlist.append(alnfile)
                tmp_files.append(alnfile)
            substrlist.append((key, alnlist))
        mican_result.append({'modeltype': modeltype, 'model': model, 'list': substrlist})
    return (mican_result, tmp_files) if return_tmpfiles==True else mican_result


def get_hits(mican_result):
    hit_dict = dict()
    for res in mican_result:
        modeltype = res['modeltype']
        model = res['model']
        substrs = res['list']
        hit_list = []
        for subname, alns in substrs:
            structure_hits = get_structure_hits(model[subname]['pdb'], alns)
            kdindex_hits = filter_by_kdindex(model[subname], structure_hits)
            hit_list = hit_list + list(kdindex_hits)
        hit_list = [(*(h['name'].split('@'))[0:2], h['kd_distance']) for h in hit_list]
        hit_dict[modeltype] = hit_list
    return hit_dict
