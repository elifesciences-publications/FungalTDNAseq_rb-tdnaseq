from __future__ import print_function
import sys
import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)
plt.style.use('seaborn-colorblind')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42



def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Plot histogram of inserts per Kb in genes')
    parser.add_argument("-p", "--pool", dest="poolfile", help="poolfile from Tnseq pipline")
    parser.add_argument("-g", "--genesGC", dest="genesGC", help="genes.GC file from tnseq pipeline (for genelengths)")
    parser.add_argument("-e", "--essential", dest="essential", help="essential genes")

    args = parser.parse_args()

    fileToOpen = args.poolfile
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            pFrame = pd.read_table(FileHandle,low_memory=False)
            FileHandle.close()
    except IOError:
        print("Could not read file:"+fileToOpen)
        sys.exit()
    pFrame['nearest_gene'].map(unicode)

    fileToOpen = args.genesGC
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            gFrame = pd.read_table(FileHandle,low_memory=False,index_col='locusId')
            FileHandle.close()
    except IOError:
        print("Could not read file:"+fileToOpen)
        sys.exit()
    gFrame.index = gFrame.index.map(unicode)

    fileToOpen = args.essential
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            eFrame = pd.read_table(FileHandle,low_memory=False,index_col='RTO4 ID')
            FileHandle.close()
    except IOError:
        print("Could not read file:"+fileToOpen)
        sys.exit()
    eFrame.index = eFrame.index.map(unicode)


    gFrame['length'] = gFrame['end'] - gFrame['begin']
    gFrame['inserts'] = 0
    gFrame['essentialOrtholog'] = 0   
    
    for geneId in gFrame.index.values:
        for geneFraction in pFrame[pFrame['nearest_gene'] == geneId]['CDS_fraction'].values:
            if 0 < int(geneFraction) < 100:
                gFrame.loc[geneId,'inserts'] += 1
        if eFrame.loc[geneId,'TotEssential'] > 0:
            gFrame.loc[geneId,'essentialOrtholog'] = 1

    gFrame['insertDensity'] = gFrame['inserts']/gFrame['length'] * 1000

    allGenesDensity = gFrame['insertDensity'].values
    orthologDensity = gFrame[gFrame['essentialOrtholog'] >= 1]['insertDensity'].values
    
    pp = PdfPages('InsertDensityHistogram.pdf')
    fig, ax = plt.subplots()
    fig.set_size_inches(3.42,2)
    bins = np.arange(0,31,0.5)
    plt.hist(allGenesDensity,bins,label='All genes')
    plt.hist(orthologDensity,bins,facecolor='C2',label='Genes with an\nessential ortholog')
    plt.xlabel('Inserts/Kb',fontsize=8,labelpad=1)
    plt.xticks(fontsize=8)
    plt.ylabel('Genes',fontsize=8,labelpad=1)
    plt.yticks(fontsize=8)
    plt.title('Histogram of insert density in genes',fontsize=10)
    legend = plt.legend(loc='upper right',edgecolor=None,frameon=False, fontsize=8)
    plt.gcf().subplots_adjust(bottom=0.16, left=0.13, right=0.98, top=0.90)
    pp.savefig()
    plt.clf()
    plt.close()
    pp.close()

    gFrame.to_csv('InsertDensitySummary.txt',sep='\t',index_label='ID')


if __name__ == "__main__":
    main(sys.argv[1:])
exit()
