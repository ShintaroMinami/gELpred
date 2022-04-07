import os

CUTOFF_RMSD = 3.0
CUTOFF_SSECOVERAGE = 0.5

def read_aln(file, name=None):
    name = os.path.splitext(os.path.basename(file))[0] if name == None else name
    with open(file, 'r') as fh:
        lines = fh.readlines()
    data, aln = [], []
    for l in lines:
        if l.startswith('# RMSD'):
            for i, rmsd in enumerate(l.split()[2:]):
                data.append(dict())
                aln.append('')
                data[i]['rmsd'] = float(rmsd)
                data[i]['sequence'] = ''
                data[i]['name'] = name+'-{:d}'.format(i+1)
        elif not l.startswith('#'):
            values = l.split()
            for i in range(len(data)):
                aln[i] = values[3*i+4]
            for i in range(len(data)):
                data[i]['sequence'] = data[i]['sequence'] + aln[i]
    return data


def read_secstr_substr(file):
    with open(file, 'r') as fh:
        lines = fh.readlines()
    for l in lines:
        if l.startswith('#') and l.endswith("SSEDEF\n"):
            ssedef = (l[1:]).split()[0]
        if l.startswith('#') and l.endswith("SSENUM\n"):
            ssenum = (l[1:]).split()[0]
    return ssedef, ssenum


def check_ssecovering(sequence, target, cutoff_coverage=CUTOFF_SSECOVERAGE):
    sses = sorted(list(set(target)))
    nlen = dict(zip(sses, [0]*len(sses)))
    nhit = nlen.copy()
    for key in nlen.keys():
        nlen[key] = target.count(key)
    for i, aa in enumerate(sequence):
        if aa != '.':
            iss = target[i]
            nhit[iss] = nhit[iss] + 1
    for key in nlen.keys():
        if key == '0': continue
        if nhit[key]/nlen[key] < cutoff_coverage:
            return False
    return True


def get_structure_hits(file_substr, list_file_aln, cutoff_rmsd=CUTOFF_RMSD):
    _, ssenum = read_secstr_substr(file_substr)
    filtered_alns = []
    for file_aln in list_file_aln:
        align = read_aln(file_aln)
        for i, aln in enumerate(align):
            i = i+1
            if aln['rmsd'] >= cutoff_rmsd: continue
            if check_ssecovering(aln['sequence'], ssenum):
                filtered_alns.append(aln)
    return sorted(filtered_alns, key=lambda x: x['name'])


def outline(file_substr, list_file_aln):
    ssedef, ssenum = read_secstr_substr(file_substr)
    alns = get_structure_hits(file_substr, list_file_aln)
    output = ""
    output += "#{}  SSES\n".format(ssedef)
    output += "#{}  SSENUM\n".format(ssenum)
    output += "#\n"
    for aln in alns:
        output += " {}  {}  {:.3f}\n".format(aln['sequence'], aln['name'], float(aln['rmsd']))
    return output
