from __future__ import print_function
import sys
import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
import re
import collections


pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)
plt.style.use('seaborn-colorblind')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42


def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Analyze biases in T-DNA insertion rates')
    parser.add_argument("-i", "--input", dest="input", help="poolfile")
    parser.add_argument("-f", "--fasta", dest="fasta", help="scaffolds in fasta format")
    parser.add_argument("-g", "--genes", dest="genes", help="position-sorted gff describing genes and rnas")

    args = parser.parse_args()

    try:
        with open(args.input, 'rb') as FileHandle:
            dFrame = pd.read_table(FileHandle,low_memory=False)
            FileHandle.close()
    except IOError:
        print("Could not read file:"+args.input)
        sys.exit()

    #Load genome
    scaffoldSequences = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    #Load gene locations
    try:
        with open(args.genes, 'rb') as FileHandle:
            genestext = FileHandle.readlines()
            FileHandle.close()
    except IOError:
        print("Could not read file:"+args.genes)
        sys.exit()

    #Order scaffolds and get lengths
    scaffoldOrder = range(len(scaffoldSequences))
    scaffoldLengths = {}
    for scaffold in scaffoldSequences:
        scaffoldOrder[int(scaffold.split("_")[1])-1] = scaffold
        scaffoldLengths[scaffold] = len(scaffoldSequences[scaffold].seq)

    #Count up inserts per scaffold
    insertionDict = {}
    Ninserts = []
    Lengths = []
    for scaffold in scaffoldOrder:
        insertionDict[scaffold] = dFrame[dFrame['scaffold'] == scaffold]['pos'].values
        Lengths.append(scaffoldLengths[scaffold])
        Ninserts.append(len(insertionDict[scaffold]))

    pp = PdfPages('AllChrom_N.pdf')

    fig, ax = plt.subplots()
    fig.text(0,0.99,'A',size=14,va='top')
    fig.set_size_inches(3.42,2)
    
    plt.scatter(Lengths,Ninserts,s=10)
    plt.xticks([0,1000000,2000000],['0','1','2'],fontsize=8)
    plt.xlabel('Scaffold length (Mb)',fontsize=8,labelpad=1)
    plt.yticks(np.arange(0,30001,10000),np.arange(0,30001,10000)/1000,fontsize=8)
    plt.ylabel('Inserts (Thousands)',fontsize=8,labelpad=1)
    plt.title('T-DNA inserts per scaffold',fontsize=10)
    plt.gcf().subplots_adjust(bottom=0.16, left=0.11, right=0.98, top=0.90)
    
    pp.savefig()
    plt.clf()
    plt.close()
    pp.close()



    #Extract flanking sequences
    print("Extracting sequences flanking insertions for GC content")
    scaffoldArr = dFrame['scaffold'].values
    posArr = dFrame['pos'].values
    flankingWindow = 50
    GCpcts = []
    scaffGC = []
    for idx, scaffold in enumerate(scaffoldArr):
        pos = int(posArr[idx])
        limits = [max(0,pos-flankingWindow),min(pos+flankingWindow,scaffoldLengths[scaffold])]
        localSeq = str(scaffoldSequences[scaffold].seq[limits[0]:limits[1]])
        Ccount = localSeq.count('C')
        Gcount = localSeq.count('G')
        GCpcts.append(100*(Gcount + Ccount)/len(localSeq))

    for scaffold in scaffoldOrder:
        randomLocations = np.random.randint(0, high=scaffoldLengths[scaffold], size=len(insertionDict[scaffold]))
        for pos in randomLocations:
            limits = [max(0,pos-flankingWindow),min(pos+flankingWindow,scaffoldLengths[scaffold])]
            localSeq = str(scaffoldSequences[scaffold].seq[limits[0]:limits[1]])
            Ccount = localSeq.count('C')
            Gcount = localSeq.count('G')
            scaffGC.append(100*(Gcount + Ccount)/len(localSeq))


    pp = PdfPages('GChistogram.pdf')

    fig, ax = plt.subplots()
    fig.text(0,0.99,'B',size=14,va='top')
    fig.set_size_inches(3.42,2)

    bins = np.arange(40,80,1)
    plt.hist(scaffGC,bins, normed=True,label='Random sites',color='C0',alpha=0.5)
    plt.hist(GCpcts,bins, normed=True,label='Insertion sites',color='C2',alpha=0.5)
    
    plt.xticks([40,50,60,70,80],fontsize=8)
    plt.xlabel('GC content (%)',fontsize=8,labelpad=1)
    plt.yticks([0,0.05,0.1],[0,5,10],fontsize=8)
    plt.ylabel('Fraction (%)',fontsize=8,labelpad=1)
    plt.title('GC content around insertion sites',fontsize=10)
    legend = plt.legend(loc='upper left',edgecolor=None,frameon=False, fontsize=8)
    plt.gcf().subplots_adjust(bottom=0.16, left=0.11, right=0.98, top=0.90)
    
    pp.savefig()
    plt.clf()
    plt.close()
    pp.close()

    


    print("Dividing genome by location type")
    #Build list of RNAs
    RNAbefore = 'telomere'
    Genes = collections.OrderedDict()
    Counter = 0
    for genesline in genestext:
        genesfeilds = genesline.split('\t')
        if len(genesfeilds) == 9:
            scaffold,evidence,featureType,start,end,flag,strand,frame,comments = genesfeilds
            start = int(start)
            end = int(end)
            if featureType == 'mRNA':
                notes = dict((k.strip(), v.strip()) for k,v in
                              (item.split('=') for item in comments.split(';')))
                promoter = []
                terminator = []
                if RNAbefore == 'telomere':
                    if strand == '+':
                        promoter =  [1,start-1]
                    else:
                        terminator =  [1,start-1]
                elif scaffold == Genes[RNAbefore]['scaffold']:
                    before = Genes[RNAbefore]['end']
                    strandBefore = Genes[RNAbefore]['strand']

                    #Check if RNAs overlap
                    if (before > start):
                        if strand == '+':
                            promoter =  [start-1,start-1]
                        else:
                            terminator =  [start-1,start-1]

                    #Otherwise, split the intergenic distance in half
                    else:
                        if strand == '+':
                            promoter =  [before+(start-before)/2,start-1]
                        else:
                            terminator =  [before+(start-before)/2,start-1]
                        if strandBefore == '+':
                            Genes[RNAbefore]['terminator'] = [before+1,before+(start-before)/2-1]
                        else:
                            Genes[RNAbefore]['promoter'] = [before+1,before+(start-before)/2-1]
                else:
                    before = Genes[RNAbefore]['end']
                    strandBefore = Genes[RNAbefore]['strand']
                    if strand == '+':
                        promoter =  [1,start-1]
                    else:
                        terminator =  [1,start-1]
                    if strandBefore == '+':
                        Genes[RNAbefore]['terminator'] = [before+1,scaffoldLengths[Genes[RNAbefore]['scaffold']]-1]
                    else:
                        Genes[RNAbefore]['promoter'] = [before+1,scaffoldLengths[Genes[RNAbefore]['scaffold']]-1]
                
                Genes[notes['ID']] = {'scaffold':scaffold,
                                             'start':int(start),
                                             'end':int(end),
                                             'strand':strand,
                                             'exons':[],
                                             'introns':[],
                                             'promoter':promoter,
                                             'terminator':terminator,
                                             '5pu':[],
                                             '3pu':[]}
                RNAbefore = notes['ID']
            elif featureType == 'CDS':
                notes = dict((k.strip(), v.strip()) for k,v in
                              (item.split('=') for item in comments.split(';')))
                if len(Genes[notes['Parent']]['exons']) > 0:
                    Genes[notes['Parent']]['introns'].append([Genes[notes['Parent']]['exons'][-1][1]+1,start-1])
                Genes[notes['Parent']]['exons'].append([start,end])
            elif featureType == 'five_prime_UTR':
                notes = dict((k.strip(), v.strip()) for k,v in
                              (item.split('=') for item in comments.split(';')))
                Genes[notes['Parent']]['5pu'] = [start,end]
            elif featureType == 'three_prime_UTR':
                notes = dict((k.strip(), v.strip()) for k,v in
                              (item.split('=') for item in comments.split(';')))
                Genes[notes['Parent']]['3pu'] = [start,end]
 
    #Complete promoter/terminator for last gene
    if Genes[RNAbefore]['strand'] == '+':
        Genes[RNAbefore]['terminator'] = [Genes[RNAbefore]['end']+1,scaffoldLengths[Genes[RNAbefore]['scaffold']]-1]
    else:
        Genes[RNAbefore]['promoter'] = [Genes[RNAbefore]['end']+1,scaffoldLengths[Genes[RNAbefore]['scaffold']]-1]


    #Correct for overlapping RNAs by splitting the overlap between UTRs
    lastGene = 'start'
    for gene in Genes:
        if (not lastGene == 'start') and (Genes[lastGene]['scaffold'] == Genes[gene]['scaffold']):
            if Genes[lastGene]['end'] > Genes[gene]['start']:
                
                    
                if Genes[gene]['strand'] == '+' and Genes[lastGene]['strand'] == '+':
                    if len(Genes[gene]['5pu']) == 0:
                        Genes[gene]['5pu'] = [Genes[gene]['start']-1,Genes[gene]['start']-1]
                    if len(Genes[lastGene]['3pu']) == 0:
                        Genes[lastGene]['3pu'] = [Genes[lastGene]['end']-1,Genes[lastGene]['end']-1]
                    UTR_start = Genes[lastGene]['3pu'][0]
                    UTR_end = Genes[gene]['5pu'][1]
                    UTR_mid = UTR_start+(UTR_start-UTR_end)/2
                    Genes[lastGene]['3pu'] = [UTR_start,UTR_mid-1]
                    Genes[gene]['5pu'] = [UTR_mid,UTR_end]
                elif Genes[gene]['strand'] == '+' and Genes[lastGene]['strand'] == '-':
                    if len(Genes[gene]['5pu']) == 0:
                        Genes[gene]['5pu'] = [Genes[gene]['start']-1,Genes[gene]['start']-1]
                    if len(Genes[lastGene]['5pu']) == 0:
                        Genes[lastGene]['5pu'] = [Genes[lastGene]['end']-1,Genes[lastGene]['end']-1]
                    UTR_start = Genes[lastGene]['5pu'][0]
                    UTR_end = Genes[gene]['5pu'][1]
                    UTR_mid = UTR_start+(UTR_start-UTR_end)/2
                    Genes[lastGene]['5pu'] = [UTR_start,UTR_mid-1]
                    Genes[gene]['5pu'] = [UTR_mid,UTR_end]
                elif Genes[gene]['strand'] == '-' and Genes[lastGene]['strand'] == '+':
                    if len(Genes[gene]['3pu']) == 0:
                        Genes[gene]['3pu'] = [Genes[gene]['start']-1,Genes[gene]['start']-1]
                    if len(Genes[lastGene]['3pu']) == 0:
                        Genes[lastGene]['3pu'] = [Genes[lastGene]['end']-1,Genes[lastGene]['end']-1]
                    UTR_start = Genes[lastGene]['3pu'][0]
                    UTR_end = Genes[gene]['3pu'][1]
                    UTR_mid = UTR_start+(UTR_start-UTR_end)/2
                    Genes[lastGene]['3pu'] = [UTR_start,UTR_mid-1]
                    Genes[gene]['3pu'] = [UTR_mid,UTR_end]
                elif Genes[gene]['strand'] == '-' and Genes[lastGene]['strand'] == '-':
                    if len(Genes[gene]['3pu']) == 0:
                        Genes[gene]['3pu'] = [Genes[gene]['start']-1,Genes[gene]['start']-1]
                    if len(Genes[lastGene]['5pu']) == 0:
                        Genes[lastGene]['5pu'] = [Genes[lastGene]['end']-1,Genes[lastGene]['end']-1]
                    UTR_start = Genes[lastGene]['5pu'][0]
                    UTR_end = Genes[gene]['3pu'][1]
                    UTR_mid = UTR_start+(UTR_start-UTR_end)/2
                    Genes[lastGene]['5pu'] = [UTR_start,UTR_mid-1]
                    Genes[gene]['3pu'] = [UTR_mid,UTR_end]
        lastGene = gene
        

    #Build list of locations in the genome and in genic and non-genic regions
    genomicLocations = {}
    promoterLocations ={}
    fivepUTRLocations = {}
    exonLocations = {}
    intronLocations = {}
    threepUTRLocations = {}
    terminatorLocations = {}
    genicLocations = {}
    intergenicLocations = {}
    for scaffold in scaffoldOrder:
        genomicLocations[scaffold] = range(1,scaffoldLengths[scaffold]+1)
        promoterLocations[scaffold] = []
        fivepUTRLocations[scaffold] = []
        exonLocations[scaffold] = []
        intronLocations[scaffold] = []
        threepUTRLocations[scaffold] = []
        terminatorLocations[scaffold] = []
        genicLocations[scaffold] = []
        intergenicLocations[scaffold] = []

    for gene in Genes:
        scaffold = Genes[gene]['scaffold']
        if len(Genes[gene]['promoter']) > 0:
            promoterLocations[scaffold].extend(range(Genes[gene]['promoter'][0],Genes[gene]['promoter'][1]+1))
        if len(Genes[gene]['terminator']) > 0:
            terminatorLocations[scaffold].extend(range(Genes[gene]['terminator'][0],Genes[gene]['terminator'][1]+1))
        for exon in Genes[gene]['exons']:
            exonLocations[scaffold].extend(range(exon[0],exon[1]+1))
        for intron in Genes[gene]['introns']:
            intronLocations[scaffold].extend(range(intron[0],intron[1]+1))
        if len(Genes[gene]['5pu']) > 0:
            fivepUTRLocations[scaffold].extend(range(Genes[gene]['5pu'][0],Genes[gene]['5pu'][1]+1))
        if len(Genes[gene]['3pu']) > 0:
            threepUTRLocations[scaffold].extend(range(Genes[gene]['3pu'][0],Genes[gene]['3pu'][1]+1))

    locationByType = {'Promoter':promoterLocations,'5pu':fivepUTRLocations,'Exon':exonLocations,'Intron':intronLocations,'3pu':threepUTRLocations,'Terminator':terminatorLocations}
    genomeDistribution = {'Promoter':0,'5pu':0,'Exon':0,'Intron':0,'3pu':0,'Terminator':0,'Total':0}
    insertionRates = {'Promoter':0,'5pu':0,'Exon':0,'Intron':0,'3pu':0,'Terminator':0,'Total':0}
    insertionClasses = ['Promoter','5pu','Exon','Intron','3pu','Terminator']
            
    for scaffold in scaffoldOrder:
        genomeDistribution['Total'] += len(genomicLocations[scaffold])
        genomeDistribution['Promoter'] += len(promoterLocations[scaffold])
        genomeDistribution['5pu'] += len(fivepUTRLocations[scaffold])
        genomeDistribution['Exon'] += len(exonLocations[scaffold])
        genomeDistribution['Intron'] += len(intronLocations[scaffold])
        genomeDistribution['3pu'] += len(threepUTRLocations[scaffold])
        genomeDistribution['Terminator'] += len(terminatorLocations[scaffold])


    #Build dictionary of insertion locations
    for index,insertInfo in dFrame.iterrows():
        if insertInfo['location'] == 'intergenic_promoter':
            insertionRates['Promoter'] += 1
        elif insertInfo['location'] == 'intergenic_terminator':
            insertionRates['Terminator'] += 1
        elif insertInfo['location'] == '5_prime_UTR':
            insertionRates['5pu'] += 1
        elif insertInfo['location'] == '3_prime_UTR':
            insertionRates['3pu'] += 1
        elif insertInfo['location'] == 'coding_region':
            if insertInfo['exon'] == 'exon':
                insertionRates['Exon'] += 1
            else:
                insertionRates['Intron'] += 1

    for insertionClass in insertionClasses:
        insertionRates['Total'] += insertionRates[insertionClass]


   





    #Simulate random
    print("Simulating biased random insertions")
    randomInsertions = {}
    biasedInsertions = {}
    for scaffold in scaffoldOrder:
        biasedInsertions[scaffold] = []
        NscaffoldInserts = len(insertionDict[scaffold])
        if NscaffoldInserts > 10 and scaffoldLengths[scaffold] > 10000:
            randomInsertions[scaffold] = np.random.randint(1, high=scaffoldLengths[scaffold]-1, size=NscaffoldInserts)
            for insertionClass in insertionClasses:
                if len(locationByType[insertionClass][scaffold]) > 10:
                    classIndexes = np.random.randint(0, high=len(locationByType[insertionClass][scaffold])-1, size=int(NscaffoldInserts*insertionRates[insertionClass]/insertionRates['Total']))
                    for i in classIndexes:
                        biasedInsertions[scaffold].append(locationByType[insertionClass][scaffold][i])

      
    

    #Taking a rolling window along first major chromosome and find the number of insertions
    window = 500
    hotspotlocations = {}
    for scaffold in scaffoldOrder[0:1]:
        hotspotlocations[scaffold] = []
        if (scaffold in biasedInsertions.keys()) and (scaffold in insertionDict.keys()):
            scaffold_length = scaffoldLengths[scaffold]

            positions = insertionDict[scaffold]

            biasedPositions = np.array(biasedInsertions[scaffold])
            locations = np.array(range(scaffold_length+window))
            scaffoldDensity = np.zeros(scaffold_length+window)
            randomDensity = np.zeros(scaffold_length+window)
            biasedDensity = np.zeros(scaffold_length+window)

            print('Summing insertions in rolling window of', 2*window, 'basepairs in ', scaffold)
            for i in locations:
                limits = [max(0,i-window),min(scaffold_length,i+window)]
                scaffoldDensity[i] = np.log10(max(((limits[0] < positions) & (positions < limits[1])).sum(),0.3))
                biasedDensity[i] = np.log10(max(((limits[0] < biasedPositions) & (biasedPositions < limits[1])).sum(),0.3))
                if scaffoldDensity[i] > 50:
                    hotspotlocations[scaffold].append(i)


            pp = PdfPages(scaffold+'_density.pdf')
            fig, ax = plt.subplots()
            fig.text(0,0.99,'C',size=14,va='top')
            fig.set_size_inches(3.42,2)
            plt.plot(locations,scaffoldDensity,label='Observed',lw=0.25)
            plt.plot(locations,biasedDensity, '-C2', label='Simulated',lw=0.25)

            if max(locations) > 1000000:
                ticksize = 500000
            elif max(locations) > 500000:
                ticksize = 100000
            elif max(locations) > 100000:
                ticksize = 50000
            elif max(locations) > 50000:
                ticksize = 10000
            else:
                ticksize = 1000
            plt.xticks(np.arange(0,max(locations),ticksize),np.arange(0,max(locations),ticksize)/1000000.0,fontsize=8)
            plt.xlabel('Position (Mb)',fontsize=8,labelpad=1)
            plt.ylabel('Inserts/Kb',fontsize=8,labelpad=1)
            plt.yticks([-0.5,0,1,2,3],[0,1,10,100],fontsize=8)
            plt.title('T-DNA insert density across scaffold '+scaffold.split("_")[1],fontsize=10)
            plt.gcf().subplots_adjust(bottom=0.16, left=0.14, right=0.98, top=0.90)
            legend = plt.legend(loc='upper center',edgecolor=None,frameon=False,fontsize=8,ncol=2)
            pp.savefig()
            plt.clf()
            plt.close()
            pp.close()



if __name__ == "__main__":
    main(sys.argv[1:])
exit()
