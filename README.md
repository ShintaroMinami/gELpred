# gELpred
GroE substrate prediction program

# Install
```bash
python setup.py install
```

# Usage
When gelpred has been installed:
```bash
gELpred  pdbfile1 pdbfile2 ...
```

W/o installation of gelpred:
```bash
./scripts/gELpred  pdbfile1 pdbfile2 ...
```

### Options
```
usage: gELpred [-h] [--model-type {RW,RWRR}] [--dir-temp DIR_TEMP]
                 pdbs [pdbs ...]

positional arguments:
  pdbs                  PDB files

optional arguments:
  -h, --help            show this help message and exit
  --model-type {RW,RWRR}, -m {RW,RWRR}
                        Prediction model type
  --dir-temp DIR_TEMP, -t DIR_TEMP
                        Temporary file directory
```