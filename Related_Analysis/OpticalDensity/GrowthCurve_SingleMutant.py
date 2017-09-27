from __future__ import print_function
import sys
import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)

def sigmoid(x, x0, k, m):
     y = m * (1 / (1 + np.exp(-k*(x-x0))))
     return y

def gsigmoid(x, K, B,v,Q):
    A=0
    #Q=0.1
    C=1
    #v=0.015
    y = A + (K - A) / ((C+Q*(np.exp(-B*x)))**(1/v))
    return y


pd.set_option("display.width",None)
pd.set_option('float_format', '{:,.2f}'.format)
plt.style.use('seaborn-colorblind')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42


def main(argv):
    parser = argparse.ArgumentParser(prog=sys.argv[0], usage='%(prog)s [options]', description='General linegraph with individual datapoints')
    parser.add_argument("-d", "--data", dest="data", help="Data table. Must have input columns Xlabel, Xvalue, Color, ColorLabel, and at least one Data column.  All columns after the first one labeled Data will be considered additional data columns")
    parser.add_argument("-o", "--output", dest="output", help="Filename for figure", default="GrowthCurve.pdf")
    parser.add_argument("-s", "--smoothing", dest="smooth", action='store_true', help="smooth averages with a sigmoidal function", default=False)
    parser.add_argument("-T", "--title", dest="title", help="title", default="Growth on Oleic Acid")
    parser.add_argument("-X", "--xlabel", dest="xlabel", help="x axis label", default="Hours")
    parser.add_argument("-Y", "--ylabel", dest="ylabel", help="y axis label", default="OD 600")
    parser.add_argument("-x", "--xticks", dest="xticks", help="x axis label", default=[0,24,48,72],nargs='+',type=int)
    parser.add_argument("-y", "--yticks", dest="yticks", help="y axis label", default=[0,2,4,6,8,10],nargs='+',type=int)
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
 

    #Start with data sorted by Xvalue
    sortedFrame = dFrame.sort_values('Xvalue')

    #Prepare figure
    figureWidth = 3.4
    figureHeight = 3
    pp = PdfPages(args.output)
    fig, ax1 = plt.subplots(figsize=(figureWidth,figureHeight))

    #Construct a non-redundant set of Xvalues and Xlabels
    Xlabeldict = {}
    Xvalues = []
    ColorLabels = []
    Colors = []
    LineStyles = []
    Symbols = []
    FaceColors = []
    for ColorLabel,colorGroup in sortedFrame.groupby('ColorOrder'):
        ColorLabels.append(colorGroup['ColorLabel'].values[0])
        Colors.append(colorGroup['Color'].values[0])
        FaceColors.append(colorGroup['FaceColor'].values[0])
        LineStyles.append(colorGroup['LineStyle'].values[0])
        Symbols.append(str(colorGroup['Symbol'].values[0]))
        for Xvalue, xGroup in colorGroup.groupby('Xvalue'):
            Xlabel = xGroup['Xlabel'].values[0]
            Xlabeldict[Xvalue] = Xlabel
            if Xvalue not in Xvalues:
                Xvalues.append(float(Xvalue))
    Xlabels = []
    for Xvalue in Xvalues:
        Xlabels.append(Xlabeldict[Xvalue])

    #Loop through colors and plot each
    Xrange = max(Xvalues) - min (Xvalues)
    for idx,ColorLabel in enumerate(ColorLabels):
        Xs = []
        Ys = []
        X_noscatter = []
        Averages = []
        AverageXs = []
        StdDevs = []
        Color = Colors[idx]
        FaceColor = FaceColors[idx]
        #Collect data into flat array and compute averages and stddevs
        for Xvalue, xGroup in sortedFrame[sortedFrame['ColorLabel'] == ColorLabel].groupby('Xvalue'):
            xData = xGroup.loc[:,'Data':].values.flatten()
            flatYs = xData[np.logical_not(np.isnan(xData))]
            stdDev = np.std(flatYs)
            average = np.mean(flatYs)
            Averages.append(average)
            StdDevs.append(stdDev)
            AverageXs.append(Xvalue)

            nYs = len(flatYs)
            xStep = 0.01*Xrange
            xDelta = -(nYs * xStep)/2
            for Y in flatYs:
                #Xs.append(Xvalue + Xrange*0.02*(0.5 - np.random.rand()))
                Xs.append(Xvalue + xDelta)
                X_noscatter.append(Xvalue)
                xDelta += xStep
            Ys.extend(flatYs)



        

        try:
            #popt, pcov = curve_fit(gsigmoid,AverageXs,Averages,maxfev=5000)
            popt, pcov = curve_fit(gsigmoid,X_noscatter,Ys,bounds=([0,0,0,0],[100,1,100,10000]))
            smoothedAverageXs = np.linspace(0, max(Xs), 50)
            smoothedAverages = gsigmoid(smoothedAverageXs, *popt)
        except OptimizeWarning:
            print("Optimize warning: couldn't fit logistic function for curve smoothing for",ColorLabel)
            smoothedAverageXs = AverageXs
            smoothedAverages = Averages
        except RuntimeError:
            print("couldn't fit logistic function for curve smoothing for",ColorLabel)
            smoothedAverageXs = AverageXs
            smoothedAverages = Averages


        
        if args.smooth:
            #ax1.plot(smoothedAverageXs,smoothedAverages,color = Color,alpha=0.9,label=ColorLabels[idx],ls=LineStyles[idx])
            ax1.plot(smoothedAverageXs,smoothedAverages,color = Color,alpha=0.9,ls=LineStyles[idx])
        else:
            #ax1.plot(AverageXs,Averages,color = Color,alpha=0.9,label=ColorLabels[idx],ls=LineStyles[idx])
            ax1.plot(AverageXs,Averages,color = Color,alpha=0.9,ls=LineStyles[idx])

        ax1.scatter(Xs,Ys,s=8,color = Color, marker=Symbols[idx],label=ColorLabels[idx],facecolor = FaceColor,lw=0.5)

        

    #Global settings for the plot
    lgnd = ax1.legend(loc='upper left',fontsize=8, ncol=1, bbox_to_anchor=(0, 1), frameon=True, labelspacing=0.1, scatterpoints=1, handletextpad=0.0)
    
    ax1.xaxis.set_ticks(args.xticks)
    #ax1.xaxis.set_ticklabels(Xlabels,rotation=90,fontsize=8)
    ax1.yaxis.set_tick_params(size=2,labelsize=8)
    ax1.yaxis.set_ticks(args.yticks)
    #ax1.yaxis.set_ticklabels(fontsize=8)
    ax1.xaxis.set_tick_params(size=2,labelsize=8)
    plt.ylabel(args.ylabel,fontsize=8,labelpad=1)
    plt.xlabel(args.xlabel,fontsize=8,labelpad=1)
    title = plt.title(args.title,fontsize=10,weight='bold')
    title.set_position([0.5,1.19])
    fig.text(0,0.99,args.letter,size=12,va='top',weight='bold')
    plt.gcf().subplots_adjust(bottom=0.1, left=0.12, right=0.95, top=0.9)

    pp.savefig()
    plt.clf()
    plt.close()
    pp.close()
    exit()
        

  

if __name__ == "__main__":
    main(sys.argv[1:])
exit()
