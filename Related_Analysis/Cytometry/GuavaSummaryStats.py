#!/usr/bin/python
from __future__ import print_function
import argparse
import sys
import numpy as np
import scipy.stats as stats
import glob
import matplotlib.pyplot as plt
import pandas as pd


def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Calulate various descriptive statistics from guava flow cytometry data')
    parser.add_argument('-p', '--path', dest='path', help="Path to directories of raw data. Results should be exported from inCyte in as list-mode data")
    parser.add_argument('-m', '--metadata', dest='metadata', help="Tab-delimited file containing meta-data listed by well.  Must contain a column corresponding to folder names for each guava run (specifiy with the -f flag) and 'Well' column")
    parser.add_argument('-f','--folderField', dest='folderField', help="Column with folder names for each entry", default="Date")
    
    args = parser.parse_args()

    if not args.path:
        print('No --path option given, trying working directory')
        path = ''
    else:
        path = args.path+'/'

    if not args.metadata:
        print('No metadata file specified. Please specifiy a metadata file with the -m flag')
        exit()
    else:
        metadata = args.metadata


    try:
        with open(args.metadata, 'rb') as FileHandle:
            mFrame = pd.read_table(FileHandle,low_memory=False,dtype='string')
            FileHandle.close()
    except IOError:
        print("Could not read file:"+args.input)
        sys.exit()


    Folders = mFrame[args.folderField].values
    Wells = mFrame['Well'].values

    BODIPY_Mean = []
    BODIPY_Geometric_Mean = []
    BODIPY_Harmonic_Mean = []
    BODIPY_Coefficient_of_Variation = []
    BODIPY_Coefficient_of_Variation_Logspace = []
    BODIPY_Skew = []
    BODIPY_Kurtosis = []
    Forward_Scatter_Mean = []
    Side_Scatter_Mean = []
    Cell_Counts = []

    for idx,well in enumerate(Wells):
        filepath = glob.glob(path + Folders[idx] + '/*' + well + '.CSV')[0]

        try:
            with open(filepath,"r") as guavaWellHandle:
                wellFrame = pd.read_table(guavaWellHandle,low_memory=False,sep=",")
                guavaWellHandle.close()
        except IOError:
            print("Could not read file:"+filepath)
            sys.exit()

        #Filter out low signal events (likely noise and cell debris)
        filteredFrame = wellFrame[ (wellFrame['GRN-HLin'] > 2) & (wellFrame['FSC-HLin'] > 50) & (wellFrame['SSC-HLin'] > 2)]

        #Collect summary stats
        BODIPY_Mean.append(filteredFrame['GRN-HLin'].mean())
        Cell_Counts.append(float(len(filteredFrame)) / ((filteredFrame['TIME'].max() / 10000.0) * 0.6))
        Forward_Scatter_Mean.append(filteredFrame['FSC-HLin'].mean())
        Side_Scatter_Mean.append(filteredFrame['SSC-HLin'].mean())
        BODIPY_Geometric_Mean.append(stats.gmean(filteredFrame['GRN-HLin'].values))
        BODIPY_Harmonic_Mean.append(stats.hmean(filteredFrame['GRN-HLin'].values))
        BODIPY_Coefficient_of_Variation.append(stats.variation(filteredFrame['GRN-HLin'].values))
        BODIPY_Coefficient_of_Variation_Logspace.append(stats.variation(filteredFrame['GRN-HLog'].values))
        BODIPY_Skew.append(stats.skew(filteredFrame['GRN-HLog'].values))
        BODIPY_Kurtosis.append(stats.kurtosis(filteredFrame['GRN-HLog'].values))


    mFrame['BODIPY'] = BODIPY_Mean
    mFrame['Cells/mL'] = Cell_Counts
    mFrame['FSC'] = Forward_Scatter_Mean
    mFrame['SSC'] = Side_Scatter_Mean
    mFrame['Geomean'] = BODIPY_Geometric_Mean
    mFrame['Harmean'] = BODIPY_Harmonic_Mean
    mFrame['CoV'] = BODIPY_Coefficient_of_Variation
    mFrame['CoVlog'] = BODIPY_Coefficient_of_Variation_Logspace
    mFrame['Skew'] = BODIPY_Skew
    mFrame['Kurtosis'] = BODIPY_Kurtosis

    mFrame.to_csv('BODIPY_sample_stats.txt',sep="\t",index=False)

  

 
if __name__ == "__main__":
    main(sys.argv[1:])
exit()




