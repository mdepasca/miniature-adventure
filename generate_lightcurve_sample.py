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


    p = subprocess.Popen("ls *.{DUMP,LIST}", shell=True, stdout=subprocess.PIPE,
            cwd=args.modelPath)
    listFiles = p.stdout.read()
    listFiles = listFiles.split('\n')
    listFiles.sort()
    listFiles.remove('')

    dumpTable = np.genfromtxt(args.modelPath+listFiles[0], skip_header=1, 
        names=True)

    dumpTable = pd.read_table(args.modelPath+listFiles[0], header=1, skiprows=0,
        sep=' ', na_values=('SN:', 'VARNAMES:')))
    
    dumpTable.sort(CID)

    listTable = np.genfromtxt(args.modelPath+listFiles[1])

    selected = np.random.choice(range(max(dumpTable.shape)), 
        size=args.numberEvents, replace=False)

    pd.to_csv(dumpTable[selected])
    np.savetxt(listTable[selected])

    

