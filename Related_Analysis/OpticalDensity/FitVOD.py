from __future__ import print_function
import sys
import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats


pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)
plt.style.use('seaborn-colorblind')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42


def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='Plot fitness vs OD data')
    parser.add_argument("-d", "--data", dest="data", help="Data table")
    parser.add_argument("-o", "--output", dest="output", help="Filename for figure", default="plot.pdf")
    parser.add_argument("-T", "--title", dest="title", help="title", default="Mutant validation in pure cultures")
    parser.add_argument("-X", "--xlabel", dest="xlabel", help="x axis label", default="Fitness score on oleic acid")
    parser.add_argument("-Y", "--ylabel", dest="ylabel", help="y axis label", default="Mutant growth (Log$_2$ vs $YKU70\Delta$)")
    parser.add_argument("-x", "--xticks", dest="xticks", help="x axis label", default=[-5,-4,-3,-2,-1,0,1,2])
    parser.add_argument("-y", "--yticks", dest="yticks", help="y axis label", default=[-5,-4,-3,-2,-1,0,1,2])
    parser.add_argument("-L", "--letter", dest="letter", help="Letter for panel", default="")

    args = parser.parse_args()

    if not args.data:
        parser.print_help()
        print("-d option required")
        sys.exit()

    fileToOpen = args.data
    try:
        with open(fileToOpen, 'rb') as FileHandle:
            dFrame = pd.read_table(FileHandle,low_memory=False)
            FileHandle.close()
    except IOError:
        print("Could not read file:"+fileToOpen)
        sys.exit()

    #Prepare figure
    figureWidth = 3.4
    figureHeight = 3.4
    pp = PdfPages(args.output)
    fig, ax1 = plt.subplots(figsize=(figureWidth,figureHeight))

  

    #Construct a non-redundant set of Xvalues and Xlabels
    Strains = []
    ShortNames = []
    Fitness = []
    ODratios = []
    FitReps = []
    ODReps = []
    FitnessAves = []
    ODAves = []
    for Strain,strainGroup in dFrame.groupby('Strain'):
        ShortName = strainGroup['Short Name'].values[0]
        OAfit = strainGroup['BarSeq Fitness'].values
        ODratio = strainGroup['OD Ratio'].values
        OAfitAverage = np.mean(OAfit)
        ODratioAverage = np.mean(ODratio)
        Strains.append(Strain)
        ShortNames.append(ShortName)
        Fitness.append(OAfitAverage)
        ODratios.append(ODratioAverage)
        FitnessAves.append([OAfitAverage,OAfitAverage,OAfitAverage])
        ODAves.append([ODratioAverage,ODratioAverage,ODratioAverage])
        FitReps.extend(OAfit)
        ODReps.extend(ODratio)

        ax1.plot(np.sort(OAfit),[ODratioAverage,ODratioAverage,ODratioAverage],lw=0.1,color='C2',zorder=1)
        ax1.plot([OAfitAverage,OAfitAverage,OAfitAverage],np.sort(ODratio),lw=0.1,color='C1',zorder=1)
        

    #ax1.scatter(AverageXs,Averages,s=10,color = Color, marker=Symbols[idx],label=ColorLabels[idx])
    ax1.scatter(Fitness,ODratios,marker='+',s=15,zorder=4,label="Average")
    ax1.scatter(FitReps,ODAves,marker='s',s=2,zorder=3,color='C2',label="Fitness Scores")
    ax1.scatter(FitnessAves,ODReps,marker='v',s=2,zorder=3,color='C1',label='Mutant Growth')


        

    #Global settings for the plot
    ax1.legend(loc='upper left',fontsize=8, frameon=True, labelspacing=0.1)

    ax1.xaxis.set_ticks(args.xticks)
    #ax1.xaxis.set_ticklabels(Xlabels,rotation=90,fontsize=8)
    ax1.yaxis.set_tick_params(size=2,labelsize=8)
    ax1.yaxis.set_ticks(args.yticks)
    #ax1.yaxis.set_ticklabels(fontsize=8)
    ax1.xaxis.set_tick_params(size=2,labelsize=8)
    ax1.axhline(y=0, color='grey',zorder=0,linewidth=0.5)
    ax1.axvline(x=0, color='grey',zorder=0,linewidth=0.5)
    plt.ylabel(args.ylabel,fontsize=8,labelpad=1)
    plt.xlabel(args.xlabel,fontsize=8,labelpad=1)
    title = plt.title(args.title,fontsize=10)
    #title.set_position([0.5,1.1])
    fig.text(0,0.99,args.letter,size=14,va='top')
    plt.gcf().subplots_adjust(bottom=0.1, left=0.12, right=0.95, top=0.9)


    #Try to pick coordinates for data labels that won't pile on top of eachother
    LabelXY = []
    for idx,X in enumerate(Fitness):
        Angle = 0
        Distance = 0.3
        TryX = X+Distance*np.cos(Angle)
        TryY = ODratios[idx]+Distance*np.sin(Angle)
        LabelXY.append([TryX,TryY])

    #Annotate datapoints
    for i, Label in enumerate(ShortNames):
        plt.annotate(Label, (LabelXY[i][0],LabelXY[i][1]),fontsize=8,zorder=4,fontname='Arial')

    


    

    

    pp.savefig()
    plt.clf()
    plt.close()
    pp.close()
    exit()
        

  

if __name__ == "__main__":
    main(sys.argv[1:])
exit()
