from __future__ import print_function
import sys
import pandas as pd
import argparse
import copy
from scipy import stats
import statsmodels.stats.multitest as smm



pd.set_option("display.width",None)
#pd.set_option('float_format', '{:,.2f}'.format)



def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='calculate enrichment of go classes')
    parser.add_argument("-l", "--list", dest="list", help="list of gene Ids")
    parser.add_argument("-g", "--genome", dest="genome", help="list of gene Ids and GO terms for the geneome")
    parser.add_argument("-t", "--ontology", dest="ontology", help="ontology descriptions")
    parser.add_argument("-o", "--outputfile", dest="output", help="name of output file")

    args = parser.parse_args()

    if not args.output:
        args.output = args.list

    #Load GO ontology
    file_to_open = args.ontology
    try:
        with open(file_to_open, 'rU') as FileHandle:
            lines = FileHandle.readlines()
            FileHandle.close()
    except IOError:
        print("Could not read file:"+file_to_open)
        sys.exit()

    goDefinitions = {}
    goNamespace = {}
    goAlternates = {}
    goID = 'null'
    name = 'null'
    namespace = 'null'
    for line in lines:
        if line[0:7] == 'id: GO:':
            goID = line[4:].rstrip()
        elif line[0:6] == 'name: ':
            name = line[6:].rstrip()
            goDefinitions[goID] = name
        elif line[0:11] == 'namespace: ':
            namespace = line[11:].rstrip()
            goNamespace[goID] = namespace
        elif line[0:7] == 'alt_id:':
            altID = line[8:].rstrip().rstrip(" ")
            goAlternates[altID] = goID

    file_to_open = args.list
    try:
        with open(file_to_open, 'rU') as FileHandle:
            lines = FileHandle.readlines()
            FileHandle.close()
    except IOError:
        print("Could not read file:"+file_to_open)
        sys.exit()

    file_to_open = args.genome
    try:
        with open(file_to_open, 'rb') as FileHandle:
            gFrame = pd.read_table(FileHandle,low_memory=False,names=['ProtID','GOID'],dtype={'ProtID':str,'GOID':str})
            FileHandle.close()
    except IOError:
        print("Could not read file:"+file_to_open)
        sys.exit()

    #replace alternate IDs with default IDs
    for checkID in gFrame['GOID'].unique():
        if checkID in goAlternates:
            for idx in gFrame[gFrame['GOID'] == checkID].index:
                gFrame.loc[idx,'GOID'] = goAlternates[checkID]

        

    #collect protein IDs for each GO term
    inputGoList = []
    inputGenesbyGOID = {}
    geneListLen = len(lines)
    for line in lines:
        line = line.rstrip()
        geneGOterms = gFrame[gFrame['ProtID'] == str(line)]['GOID'].values
        inputGoList.extend(geneGOterms)
        for GOterm in geneGOterms:
            if not GOterm in inputGenesbyGOID: 
                inputGenesbyGOID[GOterm] = [str(line)]
            else:
                inputGenesbyGOID[GOterm].append(str(line))

    countsTable = {}
    genomeLen = len(gFrame['ProtID'].unique())

    for GOID,group in gFrame.groupby('GOID'):
        inGenome = len(group)
        inList = inputGoList.count(GOID)
        if inList > 0:
            if GOID in inputGenesbyGOID:
                geneIds = inputGenesbyGOID[GOID]
            else:
                geneIds = ''
            if GOID in goDefinitions:
                description = goDefinitions[GOID]
                namespace = goNamespace[GOID]
            else:
                description = 'no definition found gene ontology file'
                namespace = 'none'
            countsTable[GOID] = {'inGenome':inGenome,
                                 'inList':inList,
                                 'GenomeLength':genomeLen,
                                 'ListLength':geneListLen,
                                 'pVal':stats.hypergeom.sf(inList,genomeLen,inGenome,geneListLen,loc=1),
                                 'geneIds':geneIds,
                                 'Description':description,
                                 'Namespace':namespace}

    statsFrame  = pd.DataFrame.from_dict(countsTable,orient='index')
    sortedFrame = statsFrame.sort_values('pVal')

   #process different namespaces seperately
    for namespace in ['molecular_function','biological_process','cellular_component']:
        namespaceFrame = sortedFrame[sortedFrame['Namespace'] == namespace].copy()

        #filtered out the less significant go terms of highly overlapping sets.
        geneIds = namespaceFrame['geneIds'].values
        filterFlag = []
        passed_geneIds = []
        for idx,genelist in enumerate(geneIds):
            filtered = False
            if (idx > 0) and len(genelist) > 0:
                for complist in passed_geneIds:
                    overlap = len(set(complist) & set(genelist))
                    if (len(genelist)*0.8 < overlap) or (len(complist)*0.7 < overlap):
                        filtered = True
            filterFlag.append(filtered)
            if not filtered:
                passed_geneIds.append(genelist)

        namespaceFrame['Filter'] = filterFlag
        filteredNSFrame = namespaceFrame[namespaceFrame['Filter'] == False].copy()

        
        #adjust pVals for multiple hypothesis
        pVals = namespaceFrame['pVal']
        adjusted=smm.multipletests(pVals,method='fdr_bh')
        namespaceFrame['Adjusted_pVal'] = adjusted[1]

        pVals = filteredNSFrame['pVal']
        adjusted=smm.multipletests(pVals,method='fdr_bh')
        filteredNSFrame['Adjusted_pVal'] = adjusted[1]

        #print to ouptput file
        namespaceFrame.reindex_axis(['Description','Namespace','GenomeLength','ListLength','inGenome','inList','pVal','Adjusted_pVal','geneIds'], axis=1).to_csv(namespace + "_" + args.output,sep='\t',index_label='GO term')
        filteredNSFrame.reindex_axis(['Description','Namespace','GenomeLength','ListLength','inGenome','inList','pVal','Adjusted_pVal','geneIds'], axis=1).to_csv(namespace + "_filtered_" + args.output,sep='\t',index_label='GO term')
        
        

    
    #Also output one list with all GOterms found
    pVals = statsFrame['pVal']
    adjusted=smm.multipletests(pVals,method='fdr_bh')

    statsFrame['Adjusted_pVal'] = adjusted[1]
    sortedFrame = statsFrame.sort_values('pVal')

    
    sortedFrame.reindex_axis(['Description','Namespace','GenomeLength','ListLength','inGenome','inList','pVal','Adjusted_pVal','geneIds'], axis=1).to_csv("all_"+args.output,sep='\t',index_label='GO term')

    quit()

    #Experimental Section to ID the most significant/largest GO term associated with each gene
    processFrame = sortedFrame[sortedFrame['Namespace'] == 'biological_process']
    print(processFrame)
    for line in lines:
        geneID = line.strip()
        largestGO = ''
        largestGOlen = 0
        Description = ''
        LowestSig = 1
        for GOID,row in processFrame.iterrows():
            genelist = inputGenesbyGOID[GOID]
            if geneID in genelist:
                #if len(genelist) > largestGOlen:
                if row['Adjusted_pVal'] < LowestSig:
                    LowestSig = row['Adjusted_pVal']
                    largestGO = GOID
                    largestGOlen = len(genelist)
                    Description = row['Description']
        print(geneID, largestGO, Description)
    

    quit()


 

    

if __name__ == "__main__":
    main(sys.argv[1:])
exit()
