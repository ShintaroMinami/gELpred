
from .utils import FILE_KDINDEX, MODEL_RW, MODEL_RWRR
import numpy as np

EPS = 0.00000001

def read_kdindex(file=FILE_KDINDEX):
    with open(file, 'r') as fh:
        lines = fh.readlines()
        kdind = dict()
        kdind['.'] = 0.0
        for l in lines:
            aa3, aa1, kd = l.split()
            kdind[aa3] = float(kd)
            kdind[aa1] = float(kd)
    return kdind

KDINDEX = read_kdindex(FILE_KDINDEX)


def read_clrep(file, id_in):
    kd_models = []
    with open(file, 'r') as fh:
        lines = fh.readlines()
        for l in lines:
            _, id, *kd = l.split()
            if int(id) == id_in:
                kd = [float(v) for v in kd]
                kd_models.append({'id': int(id), 'kd_values': kd})
    return kd_models


def filter_by_kdindex(substr, list_aln, nonredundant=True):
    kd_models = read_clrep(file=substr['clrep'], id_in=substr['id'])
    kd_cutoff = substr['kdcut']
    filtered_list = []
    for aln in list_aln:
        sequence = aln['sequence']
        kd_target = np.array([KDINDEX[aa] if aa in KDINDEX.keys() else 0.0 for aa in sequence])
        for kd_model in kd_models:
            id = kd_model['id']
            kd_values = np.array(kd_model['kd_values'])
            sum1 = np.mean(kd_target*kd_target)
            sum2 = np.mean(kd_values*kd_values)
            sumc = np.mean(kd_target*kd_values)
            kd_distance = 1 - sumc/(np.sqrt(sum1*sum2)) + EPS
            if kd_distance < kd_cutoff:
                aln['kd_distance'] = kd_distance
                aln['kd_model'] = kd_model
                filtered_list.append(aln)
    if nonredundant == False:
        return filtered_list
    nr_dict = dict()
    for aln in filtered_list:
        protname = aln['name'].split('-')[0]
        aln['protname'] = protname
        if protname in nr_dict.keys():
            if aln['kd_distance'] < nr_dict[protname]['kd_distance']:
                nr_dict[protname] = aln
        else:
            nr_dict[protname] = aln
    return nr_dict.values()
