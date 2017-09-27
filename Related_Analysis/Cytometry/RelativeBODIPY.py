from __future__ import print_function
import sys
import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


from Bio import SeqIO
import scipy.stats as stats


pd.set_option("display.width",None)
#pd.set_option('float_format', '{:,.2f}'.format)
plt.style.use('seaborn-colorblind')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42


def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Plot guava data for single deletions')
    parser.add_argument("-i", "--input", dest="input", help="Tab-delimited file of strain names and guava BODIPY averages")
    parser.add_argument("-l", "--list", dest="list", help="File with a list of strain names in order to be plotted")
    parser.add_argument("-r", "--ref", dest="ref", help="Reference strain")
    parser.add_argument("-o", "--output", dest="output", help="Outputname",default="BODIPYplot.pdf")
    parser.add_argument("-t", "--title", dest="title", help="Plot title",default="Deletion strains")
    parser.add_argument("-L", "--letter", dest="letter", help="Letter for panel", default="")
    parser.add_argument("-H", "--height", dest="height", help="plot Height", default=4)
    parser.add_argument("-W", "--width", dest="width", help="plot Width", default=2)
    parser.add_argument("-B", "--bottom", dest="bottom", help="plot Top margin (fraction)", default=0.08)
    parser.add_argument("-T", "--top", dest="top", help="plot Bottom margin (fraction)", default=0.92)
    parser.add_argument("-E", "--left", dest="left", help="plot lEft margin (fraction)", default=0.35)
    parser.add_argument("-R", "--right", dest="right", help="plot Right margin (fraction)", default=0.98)

    args = parser.parse_args()

    if not args.list:
        parser.print_help()
        print("-l option required")
        sys.exit()

    if not args.input:
        parser.print_help()
        print("-i option required")
        sys.exit()

    if not args.ref:
        parser.print_help()
        print("-r option required")
        sys.exit()

    fileToOpen = args.input
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            dFrame = pd.read_table(FileHandle,low_memory=False)
            FileHandle.close()
    except IOError:
        print("Could not read file:"+fileToOpen)
        sys.exit()

    fileToOpen = args.list
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            gFrame = pd.read_table(FileHandle,low_memory=False)
            FileHandle.close()
    except IOError:
        print("Could not read file:"+fileToOpen)
        sys.exit()

    filteredData = dFrame[dFrame['Drop'] != True]
    strainList = gFrame['Strain'].values
    strainNames = gFrame['Short Name'].values

    relativeSignal = []
    absoluteSignal = []
    dates = []
    ys = []
    absoluteAverage = []
    relativeAverage = []
    relativeMedian = []
    relativeSTD = []
    relativeSEM = []
    relativeEffectSize = []
    Ns = []
    Ttest = []
    strainIDs = []
    expectations = []
    tails = []

    


    

    
    for idx,strain in enumerate(strainList[::-1]):
        strainIDs.append(strain)
        referenceObservations = filteredData[filteredData['Strain'] == args.ref]
        refAve = np.mean(referenceObservations['BODIPY'].values)
        strainObservations = filteredData[filteredData['Strain'] == strain]

        strainAbsolute = []
        strainDates = []
        strainys = []
        strainRelative = []
        referenceRelative = []
        referenceAbsolute = []
        

        
        for index,observation in strainObservations.iterrows():
            observationDate = observation['Date']
            sameDate = referenceObservations[referenceObservations['Date'] == observationDate]
            dateAve = sameDate['BODIPY'].mean()

            for rep in sameDate['BODIPY'].values:
                referenceRelative.append(np.log2(rep/dateAve))
                referenceAbsolute.append(rep)
                
            strainAbsolute.append(observation['BODIPY'])
            strainRelative.append(np.log2(observation['BODIPY']/dateAve))
            strainDates.append(observationDate)
            strainys.append(idx+0.2-0.4*np.random.random())


        absoluteSignal.append(strainAbsolute)
        relativeSignal.append(strainRelative)
        dates.append(strainDates)
        ys.append(strainys)
        Ns.append(len(strainAbsolute))
        absoluteAverage.append(np.mean(strainAbsolute))
        relativeAverage.append(np.mean(strainRelative))
        relativeMedian.append(np.median(strainRelative))
        relativeSTD.append(np.std(strainRelative))
        relativeSEM.append(stats.sem(strainRelative))
        relativeEffectSize.append(abs(np.mean(strainRelative) - np.mean(referenceRelative))/np.sqrt(((np.std(strainRelative)**2+np.std(referenceRelative)**2)/2)))

        Tresult = stats.ttest_ind(strainRelative,referenceRelative)[1]
        expect = gFrame[gFrame['Strain'] == strain]['Expect'].values[0]


        expectations.append(expect)

        if (expect == 'none'):
            tails.append(2)
            Ttest.append(Tresult)
        else:
            tails.append(1)
            Ttest.append(Tresult/2)
            
            
    statsFrame = pd.DataFrame({'StrainID':strainIDs})
    statsFrame['N'] = Ns
    statsFrame['Mean'] = relativeAverage
    statsFrame['Median'] = relativeMedian
    statsFrame['SD'] = relativeSTD
    statsFrame['SEM'] = relativeSEM
    statsFrame['Expect'] = expectations
    statsFrame['Tails'] = tails
    statsFrame['Pval'] = Ttest
    statsFrame['D'] = relativeEffectSize
    statsFrame['Observations'] = relativeSignal

    statsFrame.to_csv(args.output+'_statsSummary.txt',sep="\t")

    flat = reduce(lambda x,y: x+y,absoluteSignal)
    flaty = reduce(lambda x,y: x+y,ys)
    flatR = reduce(lambda x,y: x+y,relativeSignal)



    pp = PdfPages(args.output+'.pdf')
    #plt.figure(figsize=(2, 5))
    fig, ax1 = plt.subplots(figsize=(float(args.width),float(args.height)))
    
    #ax2 = plt.scatter(flat,flaty,s=1,color = '#0072B2')
    #ax = plt.barh(range(0,len(absoluteAverage)),absoluteAverage,color = '#0072B2',alpha=0.3)
    ax2 = plt.plot([0,0],[-1,1000],color='grey',zorder=-3,lw=0.5)
    ax2 = plt.scatter(flatR,flaty,s=3,color = '#0072B2', alpha=0.7, edgecolor="None")
    ax = plt.barh(range(0,len(relativeAverage)),relativeAverage,color = '#0072B2',alpha=0.3)

    ax1.yaxis.set_ticks(range(len(strainNames)))
    ax1.yaxis.set_ticklabels(strainNames[::-1])
    ax1.yaxis.set_tick_params(size=2,labelsize=8)
    plt.xlabel('Relative BODIPY (log$_2$ vs WT)',fontsize=8,labelpad=1)
    plt.xticks(fontsize=8)
    ax1.set_xlim([np.nanmin(flatR)-0.5, np.nanmax(flatR)+0.5])
    ax1.set_ylim([-1,len(strainNames)])
    plt.title(args.title,fontsize=10)
    fig.text(0,0.99,args.letter,size=14,va='top')

   

    for idx,sig in enumerate(Ttest):
        if not ((relativeAverage[idx] < 0 and expectations[idx] == 'high') or (relativeAverage[idx] > 0 and expectations[idx] == 'low')):
            if sig < 0.01:
                plt.scatter(np.nanmax(flatR)+0.3,idx-0.1,marker='*',s=3,c='black')
                plt.scatter(np.nanmax(flatR)+0.3,idx+0.1,marker='*',s=3,c='black')
            elif sig < 0.05:
                plt.scatter(np.nanmax(flatR)+0.3,idx,marker='*',s=3,c='black')
    
    plt.gcf().subplots_adjust(bottom=float(args.bottom), left=float(args.left), right=float(args.right), top=float(args.top))
    pp.savefig()
    plt.clf()
    plt.close()
    pp.close()
    exit()



if __name__ == "__main__":
    main(sys.argv[1:])
exit()
