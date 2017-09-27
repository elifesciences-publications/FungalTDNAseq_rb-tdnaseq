#!/usr/bin/python
from __future__ import print_function
import sys
from optparse import OptionParser
import pandas as pd
import numpy as np
from scipy import stats
from scipy import spatial
import argparse
import statsmodels.stats.multitest as smm

np.set_printoptions(precision=2)
np.set_printoptions(linewidth=200)
pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)

def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Calculate Wilcoxon signed rank tests for barseq experiments')
    parser.add_argument("-p", "--poolcount", dest="poolCountFileName", help="PoolCount file from a TnSeq/BarSeq experiment")
    parser.add_argument("-c", "--conditions", dest="conditions", help="List of conditions and replicates formatted as a python dictionary, e.g : Condition1:[RepA,RepB,RepC],Condition2:[RepA,RepB,RepC]  Note that the order of biological replicates matters")
    parser.add_argument("-t", "--tests", dest="tests", help="Nested list of conditions to compare, eg: [[Condition1,Condition2],[Condition2,Condition3]]  Note the statistical test is only valid if replicates in the two conditions to compare are explicitly paired in some way (same preculture, paired treatment, etc).")
    parser.add_argument("-o", "--output", dest="outputbase", help="Basename for output files", default="OUT")
    parser.add_argument("-m", "--mincounts", dest="minimum_counts", help="Each insertion must have at least one condition with this many raw counts to be used.  Default is 3", default=3)
    parser.add_argument("-b", "--barcodeLocations", dest="barcodeLocationFileName", help="File with barcode locations. Columns 'nearest_gene' and 'CDS_fraction' required.  If this file is passed genewise fitness will be computed by averaging enrichment across inserts in each gene.")

    
    args = parser.parse_args()

    if not args.poolCountFileName:
        parser.print_help()
        parser.error('-p option required to select poolcount file')

    if not args.conditions:
        parser.print_help()
        parser.error('-c option required to identify conditions and replicates')

    if not args.tests:
        parser.print_help()
        parser.error('-t option required to set conditions to test')

    if not args.barcodeLocationFileName:
        parser.print_help()
        parser.error('-b option required')

    #Load poolcount file as pandas dataframe
    poolCountFileName = args.poolCountFileName
    try:
        with open(poolCountFileName, 'rb') as poolCountFileHandle:
            dataFrame = pd.read_table(poolCountFileHandle,low_memory=False)
            poolCountFileHandle.close()
    except IOError:
        print("Could not read file:"+poolCountFileName)
        sys.exit()

    #parse conditions and tests
    conditionDict = {}
    try:
        conditionString = args.conditions.replace('],',']|')
        conditionArr = conditionString.split("|")
        for conditionEntry in conditionArr:
            conditionName, repString = conditionEntry.split(":")
            repArray = repString.replace("[","").replace("]","").strip().split(",")
            conditionDict[conditionName] = repArray
    except:
        print("problem parsing -c option")

    testConditions = []
    try:
        for testString in args.tests.replace('],[','|').replace("[","").replace("]","").strip().split("|"):
            testConditions.append(testString.split(','))
    except:
        print("problem parsing -t option")

    libraryList = []
    for condition in conditionDict:
        libraryList.extend(conditionDict[condition])

        
    #Sort by location
    dataFrame.sort_values(by=['scaffold','pos'],ascending=[1,1],inplace=True)


    #Normalize by total mapped insertions per library
    print("Normalizing raw counts by total mapped insertions per library")
    normalization_factors = dataFrame.iloc[:,5:].sum().mean()/dataFrame.iloc[:,5:].sum()
    dataFrame.iloc[:,5:] = dataFrame.iloc[:,5:]*normalization_factors

    #Load barcode gene locations
    barcodeLocationFileName = args.barcodeLocationFileName
    try:
        with open(barcodeLocationFileName, 'rb') as barcodeLocationFileHandle:
            barcodeFrame = pd.read_table(barcodeLocationFileHandle,low_memory=False)
            barcodeLocationFileHandle.close()
    except IOError:
        print("Could not read file:"+barcodeLocationFileName)
        sys.exit()
    

    #Select counts data for genic barcodes and group by gene
    countsWithPositions = pd.merge(dataFrame, barcodeFrame[['barcode','nearest_gene','CDS_fraction']], how='left', on='barcode')
    central80Counts = countsWithPositions.loc[(countsWithPositions['CDS_fraction'] >= 10) & (countsWithPositions['CDS_fraction'] <= 90)]
    groupedCounts = central80Counts.groupby('nearest_gene')

    #Perfrom wilcoxn signed rank test
    print("Performing Wilcoxon Signed Rank Tests")

    geneNames = {}
    for gene, group in groupedCounts:
        nInsertions = len(group.index)
        geneNames[gene] = {'nInsertions':nInsertions}
    conditionsFrame = pd.DataFrame.from_dict(geneNames, orient='index')

    adjusted_columns = []
    for conditions in testConditions:
        comparisonLabel = ":".join(conditions)+'_wilcox'
        comparisonNranks = ":".join(conditions)+'_Nranks'
        conditionsFrame[comparisonNranks] = ''
        conditionsFrame[comparisonLabel] = ''

        for gene, group in groupedCounts:
            condition1_values = group[conditionDict[conditions[0]]].values.flatten()
            condition2_values = group[conditionDict[conditions[1]]].values.flatten()

            dropList = []
            for idx,value in enumerate(condition1_values):
                if max(condition1_values[idx],condition2_values[idx]) < float(args.minimum_counts):
                    dropList.append(idx)
            
            numRanks = np.count_nonzero(condition1_values-condition2_values)
            totalCounts = np.sum(condition1_values) + np.sum(condition2_values)
            if (numRanks > 10) and (totalCounts > 30):
                wilcox = stats.wilcoxon(condition1_values,condition2_values)
                #adust p-value below 1 to distinguish genes withoug enough inserts and reads to test
                wilcox_pval = min(wilcox[1],0.999) 
            else:
                wilcox_pval = 1
                numRanks = 0
            conditionsFrame.loc[gene,comparisonLabel] = wilcox_pval
            conditionsFrame.loc[gene,comparisonNranks] = numRanks

        #Adust for multiple hypothesis, not conting genes without enough insertions for wilcoxon test
        conditionsFrame = conditionsFrame.sort_values(by=comparisonLabel)
        uncorrected_pVals = conditionsFrame[conditionsFrame[comparisonLabel] < 1][comparisonLabel].values
        adjusted=smm.multipletests(uncorrected_pVals,method='fdr_bh')
        adjusted_pVals = np.concatenate((adjusted[1],(np.ones(len(conditionsFrame)-len(uncorrected_pVals)))),axis=0)
        conditionsFrame[comparisonLabel+'_adjusted'] = adjusted_pVals
        adjusted_columns.append(comparisonLabel+'_adjusted')

    conditionsFrame['minPval'] = conditionsFrame[adjusted_columns].min(axis=1)
    conditionsFrame = conditionsFrame.sort_values(by='minPval')

    conditionsFrame.to_csv(args.outputbase+'_Wilcoxon.txt',sep="\t",index_label='gene')

    exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
exit()




