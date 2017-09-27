#!/usr/bin/python
from __future__ import print_function
import sys
import pandas as pd
import argparse
import re

pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)

def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Combine FEBA pipeline results with effect sizes, T-like statistics and wilcoxon pvalues')
    parser.add_argument("-f", "--febadir", dest="febadir", help="Output directory from the FEBA pipeline")
    parser.add_argument("-w", "--wilcoxonfile", dest="wilcoxonfile", help="file with wilcoxon output")
    parser.add_argument("-o", "--output", dest="output", help="output filename", default="FEBA_wilcoxon.txt")
    parser.add_argument("-t", "--tests", dest="tests", help="list of conditions for which to compare T values")
    parser.add_argument("-r", "--rename", dest="rename", help="regex pattern to strip feba set prefixes from column names eg 'set[A-Z]'")

  
    args = parser.parse_args()

    if not args.febadir:
        parser.print_help()
        parser.error('FEBA output directory required. Specifiy with -f or --febadir')
    else:
        inputDirectory = args.febadir



    FEBAoutputFile = inputDirectory+'/fit_logratios.tab'
    try:
        with open(FEBAoutputFile, 'rb') as FEBAoutputFileHandle:
            fitFrame = pd.read_table(FEBAoutputFileHandle,index_col='locusId')
            FEBAoutputFileHandle.close()
    except IOError:
        print("Could not read file:"+FEBAoutputFile)
        sys.exit()

    FEBAoutputFile = inputDirectory+'/fit_t.tab'
    try:
        with open(FEBAoutputFile, 'rb') as FEBAoutputFileHandle:
            tFrame = pd.read_table(FEBAoutputFileHandle,index_col='locusId')
            FEBAoutputFileHandle.close()
    except IOError:
        print("Could not read file:"+FEBAoutputFile)
        sys.exit()


        
    tColumns = list(tFrame.columns.values)
    tColumns_renamed = tColumns[:2]
    tColumns_renamed.extend(map((lambda x: 'T_'+x),tColumns[2:]))
    tFrame.columns = tColumns_renamed

    for column in tColumns_renamed[2:]:
        fitFrame[column] = tFrame[column]

    #strip out FEBA set names
    if args.rename:
        oldColumns = list(fitFrame.columns.values)
        newColumns = oldColumns[:2]
        for column in oldColumns[2:]:
            newColumns.append(re.sub(args.rename,'',column))
        fitFrame.columns = newColumns

    #calculate new T-stats for compared conditions
    if not args.tests:
        parser.print_help()
        parser.error('-t option required to identify conditions to compare')
    else:
        testConditions = []
        try:
            for testString in args.tests.replace('],[','|').replace("[","").replace("]","").strip().split("|"):
                testConditions.append(testString.split(','))
        except:
            print("problem parsing -t option")

        for cpair in testConditions:
            T_comparison = ":".join(cpair)
            fitDiff = (fitFrame[cpair[0]] - fitFrame[cpair[1]])
            Var1 = (fitFrame[cpair[0]]/fitFrame["T_"+cpair[0]])**2.0
            Var2 = (fitFrame[cpair[1]]/fitFrame["T_"+cpair[1]])**2.0

            fitFrame[T_comparison] = fitDiff
            fitFrame["T_"+T_comparison] = fitDiff/((Var1 + Var2)**(1.0/2.0))

    #add wilcoxon signed rank test statistics
    if not args.wilcoxonfile:
        print('No wilcoxonfile supplied. Specifiy with -w or --wilcoxonfile to include wilcoxon signed rank scores in the final output')
    else:
        wilcoxonFile = args.wilcoxonfile

        FEBAoutputFile = wilcoxonFile
        try:
            with open(FEBAoutputFile, 'rb') as FEBAoutputFileHandle:
                wFrame = pd.read_table(FEBAoutputFileHandle,index_col='gene')
                FEBAoutputFileHandle.close()
        except IOError:
            print("Could not read file:"+FEBAoutputFile)
            sys.exit()

        wilcox_tests = []
        for columnName in wFrame.columns.values:
            if columnName[-15:] == 'wilcox_adjusted':
                wilcox_tests.append(columnName)

        fitFrame = fitFrame.join(wFrame[wilcox_tests])
    
    fitFrame.to_csv(args.output,sep="\t",index_label='locusId')

    exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
exit()




