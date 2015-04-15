import numpy as np
import pandas as pd
import argparse
import os
from  os import path
import subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    description = "Build sample of SNANA generated light curves, \
        comparable in numbers and type fractions to the SNPhotCC sample.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--model-path', dest='modelPath',
        help='Path where to find SNANA generated light curves using a \
            specific model (MLCS2K, SALT2, NON1A)')

    parser.add_argument(
        '-n', '--number-events', dest='numberEvents',
        default=-1, type=int,
        help='Number of SNANA simulated light curves to extract.')

    args = parser.parse_args()

else:
    pass

if __name__ == '__main__':
    if args.modelPath == '' or args.numberEvents == -1:
        raise SystemExit

    if args.modelPath[-1] != os.sep:
        args.modelPath = args.modelPath + os.sep

    finalPath = 'train_data/DES_NO_K-COR_SAMPLE/'
    if not os.path.exists(path.abspath(finalPath)):
        os.makedirs(path.abspath(finalPath))


    p = subprocess.Popen("ls *.DUMP", shell=True, stdout=subprocess.PIPE,
            cwd=args.modelPath)
    dumpFile = p.stdout.read()
    dumpFile = dumpFile.split('\n')
    dumpFile.sort()
    dumpFile.remove('')

    p = subprocess.Popen("ls *.LIST", shell=True, stdout=subprocess.PIPE,
            cwd=args.modelPath)
    listFile = p.stdout.read()
    listFile = listFile.split('\n')
    listFile.sort()
    listFile.remove('')


    dumpTable = pd.read_table(args.modelPath+dumpFile[0], header=1, skiprows=0,
        sep=' ', na_values=('SN:', 'VARNAMES:'))


    listTable = np.genfromtxt(args.modelPath+listFile[0], dtype=str)
    
    dumpTable.sort('CID')


    selected = np.random.choice(range(max(dumpTable.shape)), 
        size=args.numberEvents, replace=False)
    
    selected.sort()
    dumpTable[dumpTable.index.isin(selected)].to_csv(finalPath+dumpFile[0], 
        sep=' ', columns=('CID','GENTYPE', 'SNTYPE', 'GENZ'))

    np.savetxt(finalPath+listFile[0], listTable[selected], fmt='%s')

    

