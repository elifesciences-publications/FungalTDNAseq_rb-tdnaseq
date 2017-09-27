#!/usr/bin/python
from __future__ import print_function
import sys
import pandas as pd
import argparse

pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)

def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Average Results from FEBA barseq pipeline for biological replicates and create suitable input files for the fitness browser')
    parser.add_argument("-i", "--inputdir", dest="inputdir", help="Output directory from the FEBA pipeline")
    parser.add_argument("-o", "--outputdir", dest="outputdir", help="Directory for the output files. Defualt is to add '_averaged' to the input directory")
    parser.add_argument("-c", "--conditions", dest="conditions", help="List of conditions and replicates formatted as a python dictionary, e.g : Condition1:[RepA,RepB,RepC],Condition2:[RepA,RepB,RepC]  Note that the order of biological replicates matters")

    args = parser.parse_args()

    if not args.conditions:
        parser.print_help()
        parser.error('-c option required to identify conditions and replicates')
    else:
        conditionDict = {}
        conditionNames = []
        library_groups = []
        try:
            conditionString = args.conditions.replace('],',']|')
            conditionArr = conditionString.split("|")
            for conditionEntry in conditionArr:
                conditionName, repString = conditionEntry.split(":")
                repArray = repString.replace("[","").replace("]","").strip().split(",")
                conditionDict[conditionName] = repArray
                conditionNames.append(conditionName)
                library_groups.append(repArray)
        except:
            print("problem parsing -c option")

    if not args.inputdir:
        parser.print_help()
        parser.error('Input directory required. Specifiy input directory with -i or --inputdir')
    else:
        inputDirectory = args.inputdir

    if not args.outputdir:
        outputDirectory = inputDirectory + '_averaged'
    else:
        outputDirectory = args.outputdir

    #Pair down metafiles to reduced sample lists for averaged conditions
    filelist = {'exps':'Index',
                'expsUsed':'Index',
                'fit_quality.tab':'name'}
    for metafile in filelist:
        indexColumn = filelist[metafile]
        FEBAoutputFile = inputDirectory+'/'+metafile
        try:
            with open(FEBAoutputFile, 'rb') as FEBAoutputFileHandle:
                dataFrame = pd.read_table(FEBAoutputFileHandle)
                FEBAoutputFileHandle.close()
        except IOError:
            print("Could not read file:"+FEBAoutputFile)
            sys.exit()

        firstReps = {}
        mergelist = []
        for condition in conditionNames:
            firstReps[condition] = conditionDict[condition][0]
            mergelist.append(dataFrame[dataFrame[indexColumn] == conditionDict[condition][0]])
        
        newFrame = pd.concat(mergelist,axis=0)

        for condition in conditionNames:
            old_names = []
            rowSearch = newFrame[newFrame[indexColumn] == firstReps[condition]][indexColumn].index
            if len(rowSearch) > 0: 
                rowIndex = newFrame[newFrame[indexColumn] == firstReps[condition]][indexColumn].index.values[0]
                rowNum = newFrame.index.get_loc(rowIndex)
                columnIndex = newFrame.columns.get_loc(indexColumn)
                newFrame.iloc[rowNum,columnIndex] = condition

        if (metafile == 'expsUsed'):
            setnames = newFrame['SetName'].values
            indexnames = newFrame['Index'].values
            names = []
            for nidx, setname in enumerate(setnames):
                names.append(setname+indexnames[nidx])
            newFrame['name'] = names 
        newFrame.to_csv(outputDirectory+'/'+metafile,sep="\t",index=False)

    #Bypass single replicate quality filters and calculate the most conservative values for quality scores.
    FEBAoutputFile = inputDirectory+'/fit_quality.tab'
    try:
        with open(FEBAoutputFile, 'rb') as FEBAoutputFileHandle:
            dataFrame = pd.read_table(FEBAoutputFileHandle)
            FEBAoutputFileHandle.close()
    except IOError:
        print("Could not read file:"+FEBAoutputFile)
        sys.exit()

    setNames = dataFrame['name'].values

    setNametoRep = {}
    setNametoCond = {}
    firstSetNames = []
    conditionSetNames = {}
    for setName in setNames:
        for condition in conditionDict:
            for rep in conditionDict[condition]:
                if setName[-len(rep):] == rep:
                    setNametoRep[setName] = rep
                    setNametoCond[setName] = condition
                    if rep == conditionDict[condition][0]:
                        firstSetNames.append(setName)

    newFrame = dataFrame[dataFrame['name'].isin(firstSetNames)].copy()
        
    setConditions = {}
    for index,row in newFrame.iterrows():
        setName = row['name']
        rep = setNametoRep[setName]
        condition = setNametoCond[setName]
        newFrame.loc[index,'name'] = setName[:-len(rep)] + condition
        newFrame['u'] = 'TRUE'
        reps = conditionDict[condition]
        setNamesInCondition = []
        for setNametoCheck in setNametoRep:
            if setNametoRep[setNametoCheck] in reps:
                setNamesInCondition.append(setNametoCheck)

        setConditions[newFrame.loc[index,'name']] = setNamesInCondition
        conditionGroup = dataFrame[dataFrame['name'].isin(setNamesInCondition)]

        summedScores = ['nMapped','nPastEnd','nGenic','nUsed']
        minScores = ['gMed','gMedt0','gMean','cor12','opcor']
        maxScores = ['mad12','mad12c','mad12c_t0','adjcor','maxFit']
        maxAbsScores = ['gccor']
        for score in summedScores:
            newFrame.loc[index,score] = conditionGroup[score].sum()
        for score in minScores:
            newFrame.loc[index,score] = conditionGroup[score].min()
        for score in maxScores:
            highest = conditionGroup[score].max()
            lowest = conditionGroup[score].min()
        if -1*lowest > highest:
            newFrame.loc[index,score] = lowest
        else:
            newFrame.loc[index,score] = highest

    newFrame.to_csv(outputDirectory+'/'+metafile,sep="\t",index=False)
    conditionsToAverage = newFrame['name'].values

    #Find columns for biological replicates and average them
    filedict = {'fit_logratios.tab':{'infoColumns':['locusId','sysName','desc'],'method':1},
                'fit_logratios_unnormalized_naive.tab':{'infoColumns':['locusId','sysName','desc'],'method':1},
                'fit_logratios_unnormalized.tab':{'infoColumns':['locusId','sysName','desc'],'method':1},
                'strain_fit.tab':{'infoColumns':['barcode','rcbarcode','scaffold','strand','pos','locusId','f','used','enoughT0'],'method':1},
                'fit_t.tab':{'infoColumns':['locusId','sysName','desc'],'method':2}}
    for datafile in filedict:
        FEBAoutputFile = inputDirectory+'/'+datafile
        try:
            with open(FEBAoutputFile, 'rb') as FEBAoutputFileHandle:
                dataFrame = pd.read_table(FEBAoutputFileHandle, dtype={'locusId':str})
                FEBAoutputFileHandle.close()
        except IOError:
            print("Could not read file:"+FEBAoutputFile)
            sys.exit()

        newFrame = dataFrame[filedict[datafile]['infoColumns']].copy()
        
        columnNames = dataFrame.columns.values
        for conditionName in conditionsToAverage:
            repNames = []
            for repName in setConditions[conditionName]:
                newName = 'Match not found'
                for columnName in columnNames:
                    if columnName[:len(repName)] == repName:
                        newName = columnName
                repNames.append(newName)
            if not 'Match not found' in repNames:
                if filedict[datafile]['method'] == 1:
                    newFrame[conditionName] = dataFrame[repNames].mean(axis=1)
                elif filedict[datafile]['method'] == 2:
                    newFrame[conditionName] = dataFrame[repNames].sum(axis=1) / ((len(repNames)) ** (0.5))

        newFrame.to_csv(outputDirectory+'/'+datafile,sep="\t",index=False)

        #special case for fit_logratios.tab, save fit_logratios_good.tab bypassing single-rep quality filters
        if datafile == 'fit_logratios.tab':
            #newFrame.insert(3, "comb", concat(newFrame['sysName'],newFrame['desc']))
            newFrame.to_csv(outputDirectory+'/''fit_logratios_good.tab',sep="\t",index=False)

    exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
exit()




