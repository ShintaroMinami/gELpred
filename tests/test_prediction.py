#! /usr/bin/env python
import sys
from tqdm import tqdm
import os
maindir = os.path.abspath(os.path.dirname(__file__)+'/../')
sys.path.append(maindir)
from gepred.prediction import calc_mican, get_hits

testdir = maindir + '/tests/'

# argparse
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('model_type', type=str, choices=['RW', 'RWRR'], help='Prediction model type')
parser.add_argument('--tmp-dir', type=str, default='/var/tmp/', help='Temporary file directory')
args = parser.parse_args()

modeltype = args.model_type

file_in = testdir+'data/references/'+modeltype+'.results'
dir_test_pdb = testdir+'data/pdbs/'
dir_tmp = testdir+'tmp/'

result = dict()
with open(file_in, 'r') as fh:
    lines = fh.readlines()
    for l in lines:
        if l.startswith('#'): continue
        vals = l.split()
        name, pred = vals[0], bool(int(vals[1]))
        result[name] = {'reference': bool(int(pred))}

with tqdm(result.keys(), desc='[progress]') as pbar:
    for name in pbar:
        pdbfiles = dir_test_pdb+name+'.pdb'
        mican_result, tmp_files = calc_mican([pdbfiles], models=[modeltype], tempdir=args.tmp_dir, return_tmpfiles=True)
        # get hits
        hits = get_hits(mican_result)
        # store result
        result[name]['pred'] = True if len(hits[modeltype]) > 0 else False
        check = 'OK' if result[name]['pred'] == result[name]['reference'] else 'FAIL'
        pbar.set_postfix({name: check})

count_ok = 0
count_tot = 0
for k in result.keys():
    count_tot += 1
    p0 = result[k]['reference']
    p = result[k]['pred']
    if(p0 == p):
        count_ok += 1
        continue
    else:
        print('ERROR: ', k, p0, p)
print(" SUCCESS: {} / {}".format(count_ok, count_tot))

for f in tmp_files:
    os.remove(f)