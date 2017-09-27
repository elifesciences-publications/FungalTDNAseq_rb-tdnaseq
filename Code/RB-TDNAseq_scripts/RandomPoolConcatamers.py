#!/usr/bin/python
from __future__ import print_function
import sys
from optparse import OptionParser
import numpy as np
import operator
import pandas as pd
from datetime import datetime

#http://code.activestate.com/recipes/576874-levenshtein-distance/
def levenshtein(s1, s2):
    l1 = len(s1)
    l2 = len(s2)
    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in range(l2 + 1):
        matrix[zz] = range(zz,zz + l1 + 1)
    for zz in range(0,l2):
        for sz in range(0,l1):
            if s1[sz] == s2[zz]:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
            else:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]

def ReverseComplement(seq):
    # too lazy to construct the dictionary manually, use a dict comprehension
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = {}
    for i in range(16):
        if i < 4 or 8<=i<12:
            seq_dict[seq1[i]] = seq1[i+4]
    #seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])


def OffByOneList(seq):
    if seq[0] in ("A","T","G","C"):
        char_set = ("A","T","G","C")
    elif seq[0] in ("a","t","g","c"):
        char_set = ("a","t","g","c")
    else:
        return False

    variants = []
    for chari in range(len(seq)):
        if chari == 0:
            preseq = ""
        else:
            preseq = seq[0:chari]
        if chari == len(seq)-1:
            postseq = ""
        else:
            postseq = seq[chari+1:]
        for char in char_set:
            if seq[chari] != char:
                variants.append(preseq+char+postseq)

    return(variants)

def main(argv):
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="mappedFilesList", help="Comma seperated list of mapped reads (TAB files from TnSeq mapping step)")
    parser.add_option("-f", "--minFrac", dest="minFraction", help="Minimum fraction of reads mapping to a single location", default=0.8)
    parser.add_option("-m", "--maxQbeg", dest="maxQbeg", help="Maximum number of bases after junction before start of match to genome", default=3)
    parser.add_option("-o", "--output", dest="outPutFileName", help="Basename for output files", default="pool")
    parser.add_option("-t", "--offbytwo", dest="offbytwo", action='store_true', help="optional filter of off-by-two barcodes", default=False)
    (options, args) = parser.parse_args()

    if not options.mappedFilesList:
        parser.print_help()
        parser.error('Input file required. Specifiy input file with -i or --input')
        


    totals={}
    barcode_to_genome = {}
    barcode_to_pastEnd = {}
    barcode_flags = {}
    barcode_primary_location = {}
    barcode_summaries = {}
    offbytwo = []

    #Load up mapped inserts from one or more files
    tabFileNames = options.mappedFilesList
    minFraction = float(options.minFraction)
    maxQbeg = int(options.maxQbeg)


    totalLines = 0
    nFiles = 0
    for tabFileName in tabFileNames.split(","):
        nLines = 0
        timestamp = "["+datetime.now().time().strftime('%X')+"] "
        print(timestamp+"Reading tabfile: "+tabFileName)
        try:
            with open(tabFileName, 'rb') as tabFileHandle:

                for line in tabFileHandle:
                    nLines += 1
                    if (nLines % 1000000 == 0 ):
                        timestamp = "["+datetime.now().time().strftime('%X')+"] "
                        print(timestamp, (nLines/1000000), " million lines read...")

                    readMapped = line.rstrip().split()
                    location = readMapped[2]+":"+readMapped[3]+":"+readMapped[4]

                    #First check if barcode is new or has been seen before
                    newMappedBarcode = readMapped[1]
                    if "N" not in newMappedBarcode:
                        if newMappedBarcode in totals:
                                totals[newMappedBarcode] += 1
                                if readMapped[2] == "pastEnd":
                                    barcode_flags[newMappedBarcode]['pastEndTot'] += 1
                                    if location in barcode_to_pastEnd[newMappedBarcode]:
                                        barcode_to_pastEnd[newMappedBarcode][location] += 1
                                    else:
                                        barcode_to_pastEnd[newMappedBarcode][location] = 1
                                else:
                                    if location in barcode_to_genome[newMappedBarcode]:
                                        barcode_to_genome[newMappedBarcode][location] += 1
                                    else:
                                        barcode_to_genome[newMappedBarcode][location] = 1
                        else:
                                totals[newMappedBarcode] = 1
                                if readMapped[2] == "pastEnd":
                                    barcode_to_pastEnd[newMappedBarcode] = {location:1}
                                    barcode_to_genome[newMappedBarcode] = {'Null:0:+':0}
                                    barcode_flags[newMappedBarcode] = {'ambReads':0,'highQbeg':0,'pastEndTot':1}
                                else:
                                    barcode_to_genome[newMappedBarcode] = {location:1}
                                    barcode_to_pastEnd[newMappedBarcode] = {'Null:0:+':0}
                                    barcode_flags[newMappedBarcode] = {'ambReads':0,'highQbeg':0,'pastEndTot':0}

                        if int(readMapped[6]) > maxQbeg:
                            barcode_flags[newMappedBarcode]['highQbeg'] += 1

                        if (int(readMapped[5]) != 1):
                            barcode_flags[newMappedBarcode]['ambReads'] += 1

                tabFileHandle.close()
                totalLines = totalLines + nLines
                nFiles += 1

        except IOError:
            print("Could not read file:"+tabFileName)
            sys.exit()

    print(str(totalLines)+" lines read from "+str(nFiles)+" files")


    # Mask out off-by-one variants
    timestamp = "["+datetime.now().time().strftime('%X')+"] "
    print(timestamp+"Masking out off-by-one barcodes. (Most likely sequencing errors)")

    offByOneList = []
    for barcode in totals:
        variants = OffByOneList(barcode)
        offByOne = False
        for variantBarcode in variants:
            if (not offByOne) & (variantBarcode in totals):
                if (totals[variantBarcode] > totals[barcode]*100):
                    offByOne = True
                    offByOneList.append(barcode)

    outPutFileName = options.outPutFileName
    offByOneFileHandle = open(outPutFileName+"_offByOne", 'w')
    barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n"])
    offByOneFileHandle.write(barcodeSummary+"\n")
    nOffByOne = 0
    OffByOneReads = 0
    for barcode in offByOneList:
        OffByOneReads += totals[barcode]
        barcodeSummary = "\t".join([barcode,ReverseComplement(barcode),str(totals[barcode]), str(totals[barcode]-barcode_flags[barcode]['highQbeg'])])
        offByOneFileHandle.write(barcodeSummary+"\n")
        del barcode_to_pastEnd[barcode]
        del barcode_to_genome[barcode]
        del totals[barcode]
        del barcode_flags[barcode]
        nOffByOne +=1
    offByOneFileHandle.close()

    timestamp = "["+datetime.now().time().strftime('%X')+"] "
    print(timestamp+"Wrote", nOffByOne, "off-by-one barcodes ("+str(OffByOneReads)+" reads) to "+outPutFileName+"_offByOne")


    #Set up pool file
    outPutFileName = options.outPutFileName
    outPutFileHandle = open(outPutFileName, 'w')
    ambFileHandle = open(outPutFileName+"_ambiguous", 'w')
    multiFileHandle = open(outPutFileName+"_multiLocus", 'w')
    barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n", "scaffold", "strand", "pos","type","nMainLocation", "nTot", "nPastEnd", "All genomic mappings", "All pastEnd mappings"])
    outPutFileHandle.write(barcodeSummary+"\n")
    ambFileHandle.write(barcodeSummary+"\n")
    multiFileHandle.write(barcodeSummary+"\n")


    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp,"Processing barcodes with mapped genomic locations")

    nMulti = 0
    nAmb = 0
    nPastEnd = 0
    nMapped = 0
    MultiReads = 0
    AmbReads = 0
    allPastEndReads = 0
    MappedReads = 0
    
    for barcode, barcode_locations in barcode_to_genome.items():

        sorted_locations = sorted(barcode_to_genome[barcode].items(), key=operator.itemgetter(1), reverse=True)
        firstLocation = sorted_locations[0][0].split(":")

        ambiguousReads = int(barcode_flags[barcode]['ambReads'])
        pastEndReads = int(barcode_flags[barcode]['pastEndTot'])
        highQbegreads = int(barcode_flags[barcode]['highQbeg'])
        insertionType = ""
        primeLocationAbundance = 0
        
        #Ignore barcodes with vast majority of reads pastEnd
        if (pastEndReads > (sorted_locations[0][1] * 7)):
            insertionType = "pastEnd"
            primeLocationAbundance = pastEndReads

        #If ambiguous reads are more common than top unambiguous location tag, barcode with ambiguous mapping
        elif (sorted_locations[0][1] < 2 * ambiguousReads):
            insertionType = "Amb"
            primeLocationAbundance = ambiguousReads

        #If not pastEnd, ambiguous, or multilocus hit, classify the type of insertion
        else:
            #If second most common location is on opposite strand nearby, call as inverted concatamer
            if (len(sorted_locations) > 1) and (sorted_locations[1][0] != 'Null'):
                secondLocation = sorted_locations[1][0].split(":")
                if (firstLocation[0] == secondLocation[0]) & (firstLocation[2] != secondLocation[2]) & (abs(int(firstLocation[1])-int(secondLocation[1])) < 10):
                    insertionType = "InCat"
                    primeLocationAbundance = sorted_locations[0][1]+sorted_locations[1][1]
                else:
                    insertionType = "Single"
                    primeLocationAbundance = sorted_locations[0][1]
            else:
                insertionType = "Single"
                primeLocationAbundance = sorted_locations[0][1]

            #If not inverted concatamer, but both pastEnd hits and genome hits found, tag as concatamer
            if (pastEndReads > (sorted_locations[0][1])/10) & (insertionType != "InCat"):
                insertionType = "Concat"
                primeLocationAbundance = sorted_locations[0][1]

            #Most common location must accout for at least minFrac of reads (not counting pastEnd)
            if ((float(primeLocationAbundance)/float(totals[barcode]-pastEndReads)) < minFraction):
                insertionType = "Multi"
                primeLocationAbundance = sorted_locations[0][1]

        #Print out summary
        barcodeSummary = "\t".join([barcode,ReverseComplement(barcode),str(totals[barcode]), str(totals[barcode]-highQbegreads), firstLocation[0], firstLocation[2], firstLocation[1],insertionType,str(primeLocationAbundance), str(totals[barcode]), str(pastEndReads), str(sorted_locations), str(barcode_to_pastEnd[barcode])])

        
        
        if insertionType in ["Single","Concat","InCat"]:
            outPutFileHandle.write(barcodeSummary+"\n")
            nMapped += 1
            MappedReads += totals[barcode]
            barcode_primary_location[barcode] = {'scaffold':firstLocation[0],'strand':firstLocation[2],'position':int(firstLocation[1]),'counts':int(totals[barcode])}
            barcode_summaries[barcode] = barcodeSummary
        elif insertionType == "Multi":
            multiFileHandle.write(barcodeSummary+"\n")
            nMulti += 1
            MultiReads += totals[barcode]
        elif insertionType == "Amb":
            ambFileHandle.write(barcodeSummary+"\n")
            nAmb += 1
            AmbReads += totals[barcode]

        #Remove barcode from pastEnd list to prevent redundant reporting
        if insertionType != "pastEnd":
            del barcode_to_pastEnd[barcode]

    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp,"Wrote",nMapped,"barcodes ("+str(MappedReads)+" reads) to",outPutFileName)
    outPutFileHandle.close()
    print(timestamp,"Wrote",nMulti,"barcodes ("+str(MultiReads)+" reads) to",outPutFileName+"_multiLocus")
    multiFileHandle.close()
    print(timestamp,"Wrote",nAmb,"barcodes ("+str(AmbReads)+" reads) to",outPutFileName+"_ambiguous")
    ambFileHandle.close()

    #Set up pool file for pastEnd hits
    outPutFileHandle = open(outPutFileName+"_pastEnd", 'w')
    barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n", "scaffold", "strand", "pos","type","nHeadtoTail", "nHeadtoHead", "nUnprocessed", "Misprocessed"])
    outPutFileHandle.write(barcodeSummary+"\n")

    timestamp = "["+datetime.now().time().strftime('%X')+"] "
    print(timestamp+"Processing pastEnd hits")

    #Report barcodes with only pastEnd reads
    for barcode in barcode_to_pastEnd:
        HeadtoTail = 0
        HeadtoHead = 0
        Unprocessed = 0
        Misprocessed = 0

        sorted_locations = sorted(barcode_to_pastEnd[barcode].items(), key=operator.itemgetter(1), reverse=True)
        firstLocation = sorted_locations[0][0].split(":")
        pastEndReads = barcode_flags[barcode]['pastEndTot']

        #Classify the type of insertion
        for pastEndLocation in sorted_locations:
            pastEndLocationDetails = pastEndLocation[0].split(":")

            if (pastEndLocationDetails[2] == '-') & (int(pastEndLocationDetails[1]) > 7000) :
                HeadtoHead += pastEndLocation[1]
            elif (pastEndLocationDetails[2] == '+') & (int(pastEndLocationDetails[1]) > 5000):
                HeadtoTail += pastEndLocation[1]
            elif (pastEndLocationDetails[2] == '+') &  (int(pastEndLocationDetails[1]) < 10):
                Unprocessed += pastEndLocation[1]
            else:
                Misprocessed += pastEndLocation[1]

        if HeadtoHead > pastEndReads * 0.9:
            insertionType = "HeadtoHead"
        elif Unprocessed > pastEndReads * 0.9 :
            insertionType = "Unprocessed"
        elif HeadtoTail > pastEndReads * 0.9:
            insertionType = "HeadtoTail"
        elif Misprocessed > pastEndReads * 0.9:
            insertionType = "Misprocessed"
        else:
            insertionType = "Mixed"


        #Print out summary
        barcodeSummary = "\t".join([barcode,ReverseComplement(barcode),str(pastEndReads), str(pastEndReads), "pastEnd", firstLocation[2], firstLocation[1],insertionType,str(HeadtoTail), str(HeadtoHead), str(Unprocessed),str(Misprocessed), str(sorted_locations), str(barcode_to_genome[barcode])])
        outPutFileHandle.write(barcodeSummary+"\n")
        nPastEnd += 1
        allPastEndReads += totals[barcode]
    outPutFileHandle.close()

    timestamp = "["+datetime.now().time().strftime('%X')+"]"
    print(timestamp,"Wrote",nPastEnd,"barcodes ("+str(allPastEndReads)+" reads) to ",outPutFileName+"_pastEnd")

   
    if options.offbytwo:
        #Find similair sequences (off-by-two, short indels) at same position (computationaly intensive).
        timestamp = "["+datetime.now().time().strftime('%X')+"]"
        print(timestamp,"Finding similair barcodes that mapped to the same location (off-by-two, indels, etc).  This is a computationally intensive step and may take a while...")

        #Loop through scaffolds
        locationFrame = pd.DataFrame.from_dict(barcode_primary_location,orient='index')
        scaffoldGroups = locationFrame.groupby('scaffold')
        NoffByTwo = 0
        for scaffold, barcodes in scaffoldGroups:
            timestamp = "["+datetime.now().time().strftime('%X')+"]"
            print(timestamp,"Checking scaffold ",scaffold)
            #For each scaffold get a list of barcodes sorted by position
            barcodes = barcodes.sort_values(by='position')
            scaffold_barcode_list = barcodes.index.values
            scaffold_barcode_positions = barcodes['position'].values
            scaffold_barcode_counts = barcodes['counts'].values

            #Scroll a moving window through each scaffold checking local neighborhood for similair barcodes
            for idx, ref_position in enumerate(scaffold_barcode_positions[0:-1]):
                window_end = min(len(scaffold_barcode_positions),idx+20)
                check_idx = idx+1
                #only need to compare barcodes that map very close to each other
                while check_idx < len(scaffold_barcode_positions) and scaffold_barcode_positions[check_idx] - ref_position < 11:
                    #can only determine which barcode is correct if one is much more abundant than the other
                    if abs(np.log2(float(scaffold_barcode_counts[idx])/float(scaffold_barcode_counts[check_idx]))) > 1:
                        #Calculate minimum edit distance for nearby by barcodes when one is much more abundant than the other
                        if levenshtein(scaffold_barcode_list[idx], scaffold_barcode_list[check_idx]) < 5:
                            NoffByTwo += 1
                            if scaffold_barcode_counts[idx] > scaffold_barcode_counts[check_idx]:
                                offbytwo.append(scaffold_barcode_list[check_idx])
                            else:
                                offbytwo.append(scaffold_barcode_list[idx])
                    check_idx += 1

        outPutFileHandle = open(outPutFileName, 'w')
        obtFileHandle = open(outPutFileName+"_offbytwo", 'w')
        barcodeSummary = "\t".join(["barcode","rcbarcode","nTot", "n", "scaffold", "strand", "pos","type","nMainLocation", "nTot", "nPastEnd", "All genomic mappings", "All pastEnd mappings"])
        obtFileHandle.write(barcodeSummary+"\n")
        outPutFileHandle.write(barcodeSummary+"\n")
        for barcode in barcode_summaries:
            summary_arr = barcode_summaries[barcode].split("\t")
            if barcode in offbytwo:
                summary_arr[7] = "offByTwo"
                obtFileHandle.write("\t".join(summary_arr)+"\n")
            else:
                outPutFileHandle.write("\t".join(summary_arr)+"\n")
        outPutFileHandle.close()
        obtFileHandle.close()
        timestamp = "["+datetime.now().time().strftime('%X')+"]"
        print(timestamp,"Moved",NoffByTwo,"barcodes with likely sequencing errors to ",outPutFileName+"_offbytwo")


if __name__ == "__main__":
    main(sys.argv[1:])
exit()




