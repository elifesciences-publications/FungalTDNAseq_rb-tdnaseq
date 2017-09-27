#!/usr/bin/python
from __future__ import print_function
import sys
from optparse import OptionParser
import pandas as pd
import numpy as np

from datetime import datetime

pd.set_option("display.width",None)

def usage():
    usage_text = "Catagorize the location of insertions in a mutant pool as intergenic,utr, intronic or exonic"
    usage_text += "\nExample: PlotRegion.py -e enrichments -g GFF_file -c scaffold:coord1:coord2 -l libraries"
    print(usage_text)

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--pool_file", dest="pool_file", help="Poolfile produced by RandomPoolConcatamers.py")
    parser.add_option("-g", "--gff", dest="gff", help="GFF3 with feature locations")
    parser.add_option("-o", "--output", dest="outputbase", help="Basename for output files", default="pool")
    #parser.add_option("-i", "--identifier", dest="identifier", help="gene identifier name, eg. for fields annotated as 'orf_end=7629; gene_name gene_1.2', the identifier could be 'orf_end=' or 'gene_name'", default='gene_name')

    (options, args) = parser.parse_args()


    if not options.pool_file or not options.gff:
        parser.print_help()
        parser.error('A pool file and gff are required')
        

    #Load enrichment file and GFF as pandas dataframes
    poolFileName = options.pool_file
    gffFileName = options.gff

    try:
        with open(gffFileName, 'rb') as gffFileHandle:
            gffContent = gffFileHandle.readlines()
            

            gffFileHandle.close()
            timestamp = "["+datetime.now().time().strftime('%X')+"]"
            print(timestamp, "Read", len( gffContent), "features from",gffFileName)
    except IOError:
        print("Could not read file:"+gffFileName)
        sys.exit()
    
    try:
        with open(poolFileName, 'rb') as poolFileHandle:
            poolFrame = pd.read_table(poolFileHandle)
            poolFileHandle.close()
            poolFrame.dropna(how='all')
            timestamp = "["+datetime.now().time().strftime('%X')+"]"
            print(timestamp, "Read", len(poolFrame), "barcodes from",poolFileName)
    except IOError:
        print("Could not read file:"+poolFileName)
        sys.exit()


    #Sort pool by location
    poolFrame.sort_values(['scaffold','pos'],ascending=[1,1],inplace=1)
    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp, "Sorted pool entries by location")
    
    #Extract arrays of scaffolds and positions to check against GFF
    scaffolds = poolFrame['scaffold'].values
    positions = poolFrame['pos'].values

    #Extract gene identifiers from GFF
    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp, "Parsing gene info from GFF")

    gffcolumns = ['scaffold','type','feature','begin','end','flag','strand','frame','attributes']
    mRNA_Dict = {}

    #First create a dictionary for every mRNA
    for gffline in gffContent:
        if not gffline.strip()[0] == "#":
            lineDict = dict(zip(gffcolumns,gffline.split("\t")))
            if lineDict['feature'] == "mRNA":
                attributes = dict(item.split("=") for item in lineDict['attributes'].split(";"))
                mRNA_Dict[attributes['ID']] = {'gene_id':attributes['proteinId'],
                                                       'mRNA':attributes['ID'],
                                                       'scaffold':lineDict['scaffold'],
                                                       'strand':('+' == lineDict['strand']),
                                                       'mRNA_begin':int(lineDict['begin']),
                                                       'mRNA_end':int(lineDict['end']),
                                                       'CDS_begin':-1,
                                                       'CDS_end':-1,
                                                       'exon_begin':[],
                                                       'exon_end':[]}
                
    #Populate the coding exons for this mRNA (not including utrs)
    for gffline in gffContent:
        if not gffline.strip()[0] == "#":
            lineDict = dict(zip(gffcolumns,gffline.split("\t")))
            if lineDict['feature'] == "CDS":
                attributes = dict(item.split("=") for item in lineDict['attributes'].split(";"))
                parent = attributes['Parent'].strip()
                mRNA_Dict[parent]['exon_begin'].append(int(lineDict['begin']))
                mRNA_Dict[parent]['exon_end'].append(int(lineDict['end']))


    #Translate to a pandas dataframe and sort by location. DataFrame constructor is choking on exon_begins/ends arrays, so need to build clumsy lists first and add as columns.
    protIDs = []
    gene_scaffolds = []
    strands = []
    mRNAs = []
    mRNA_begins = []
    mRNA_ends = []
    CDS_begins = []
    CDS_ends = []
    exon_begins = []
    exon_ends = []
    
    for mRNA in mRNA_Dict:
        protIDs.append(mRNA_Dict[mRNA]['gene_id'])
        gene_scaffolds.append(mRNA_Dict[mRNA]['scaffold'])
        strands.append(mRNA_Dict[mRNA]['strand'])
        mRNAs.append(mRNA)
        mRNA_begins.append(mRNA_Dict[mRNA]['mRNA_begin'])
        mRNA_ends.append(mRNA_Dict[mRNA]['mRNA_end'])
        CDS_begins.append(min(mRNA_Dict[mRNA]['exon_begin']))
        CDS_ends.append(max(mRNA_Dict[mRNA]['exon_end']))
        exon_begins.append(mRNA_Dict[mRNA]['exon_begin'])
        exon_ends.append(mRNA_Dict[mRNA]['exon_end'])


    geneFrame = pd.DataFrame({'gene_id':protIDs,
                              'scaffold':gene_scaffolds,
                              'strand':strands,
                              'mRNA':mRNAs,
                              'mRNA_begin':mRNA_begins,
                              'mRNA_end':mRNA_ends,
                              'CDS_begin':CDS_begins,
                              'CDS_end':CDS_ends,
                              'exon_begin':exon_begins,
                              'exon_end':exon_ends},
                              index=protIDs)


    geneFrame.sort_values(['scaffold','CDS_begin'],ascending=[1,1],inplace=1)
    geneFrame.to_csv('geneFrame.txt',sep="\t")
        
    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp, len(geneFrame), "Genes Processed")
    print(timestamp, "Processing insertions in mutant pool. (This may take several minutes)")
    
    #Loop through positions identifying features in GFF
    genic_types = []
    exonic_types = []
    nearest_genes = []
    coding_fractions = []
    bp_aways = []
    for idx,scaffold in enumerate(scaffolds):
        if idx in [1000,10000,100000]:
            timestamp = "["+datetime.now().time().strftime('%X')+"]"
            print(timestamp, idx, "Insertions Processed")

        position = positions[idx]

        scaffoldFrame = geneFrame[geneFrame['scaffold'] == scaffold]
        geneAfter = scaffoldFrame[scaffoldFrame['CDS_end'] >= position][:1]
        geneBefore = scaffoldFrame[scaffoldFrame['CDS_begin'] <= position][-1:]

        neighborhood = pd.concat([geneAfter,geneBefore])

        if len(neighborhood) == 0:
            nearest_genes.append("empty_scaffold")
            genic_types.append("intergenic")
            coding_fractions.append(-1)
            exonic_types.append("intergenic")
            bp_aways.append("NA")
            
            
        else:
            neighborhood['end_dist'] = abs(position - neighborhood['CDS_end'])
            neighborhood['begin_dist'] = abs(position - neighborhood['CDS_begin'])
            neighborhood['min_dist'] = np.where(neighborhood['end_dist'] > neighborhood['begin_dist'],neighborhood['begin_dist'],neighborhood['end_dist'])
            min_idx = np.where(neighborhood['end_dist'] > neighborhood['begin_dist'],neighborhood['begin_dist'],neighborhood['end_dist']).argmin()
            local_gene = neighborhood.iloc[min_idx].to_dict()

            nearest_gene = local_gene['gene_id']

            if local_gene['strand']:
                local_coordinates = {'TSS':local_gene['mRNA_begin']-local_gene['CDS_begin'],
                                     'Stop':local_gene['CDS_end']-local_gene['CDS_begin'],
                                     'Term':local_gene['mRNA_end']-local_gene['CDS_begin'],
                                     'Ins':position-local_gene['CDS_begin']}
            else:
                local_coordinates = {'TSS':local_gene['CDS_end']-local_gene['mRNA_end'],
                                     'Stop':local_gene['CDS_end']-local_gene['CDS_begin'],
                                     'Term':local_gene['CDS_end']-local_gene['mRNA_begin'],
                                     'Ins':local_gene['CDS_end']-position}

            if local_coordinates['Ins'] < local_coordinates['TSS']:
                genic_type = 'intergenic_promoter'
                exonic_type = 'intergenic'
                coding_fraction = -1
                bp_away = local_coordinates['Ins']
                
            elif local_coordinates['Ins'] < 0:
                genic_type = '5_prime_UTR'
                exonic_type = '5_prime_UTR'
                coding_fraction = -1
                bp_away = local_coordinates['Ins']
                
            elif local_coordinates['Ins'] <= local_coordinates['Stop']:
                genic_type = 'coding_region'
                exonic_type = 'intron'
                coding_fraction = local_coordinates['Ins'] * 100 / local_coordinates['Stop']
                bp_away = 0
                exons = zip(local_gene['exon_begin'],local_gene['exon_end'])
                for exon in exons:
                    if exon[0] <= position <= exon[1]:
                        exonic_type = 'exon'
            elif local_coordinates['Ins'] <= local_coordinates['Term']:
                genic_type = '3_prime_UTR'
                exonic_type = '3_prime_UTR'
                coding_fraction = -1
                bp_away = local_coordinates['Ins']
            else:
                genic_type = 'intergenic_terminator'
                exonic_type = 'intergenic'
                coding_fraction = -1
                bp_away = local_coordinates['Ins']

            nearest_genes.append(nearest_gene)
            genic_types.append(genic_type)
            coding_fractions.append(coding_fraction)
            exonic_types.append(exonic_type)
            bp_aways.append(bp_away)
        
    poolFrame['nearest_gene'] = nearest_genes
    poolFrame['location'] = genic_types
    poolFrame['CDS_fraction'] = coding_fractions
    poolFrame['Basepair_from_CDS'] = bp_aways
    poolFrame['exon'] = exonic_types

    poolFrame.to_csv(options.outputbase+'_genes.txt',sep="\t",index=False)
    

    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp, "Wrote poolfile with gene information to", options.outputbase+'_genes.txt')
 

if __name__ == "__main__":
    main(sys.argv[1:])
exit()




