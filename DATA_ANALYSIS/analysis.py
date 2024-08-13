import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from scipy.interpolate import make_interp_spline 
import numpy as np


def average(list):
    return sum(list)/len(list)

def standarddeviation(list):
    averageofthelist = average(list)
    newlist = [ (elementofthelist - averageofthelist)**2 for elementofthelist in list ]
    sigmasquared  = average(newlist)
    sigma = sigmasquared**0.5
    return sigma

def autocorrelation(list, lag):
    xi = list[:-lag]
    xiplust = list[lag:]
    xitimesxiplust = [ i*iplust for (i, iplust) in zip(xi, xiplust) ]
    averageofxitimexiplust = average(xitimesxiplust)
    averageofxi = average(xi)
    averageofxiplust = average(xiplust)
    autocorrelationvalue = averageofxitimexiplust - (averageofxi*averageofxiplust)
    return autocorrelationvalue

def integratedautocorrelation(list):
    selfcorrelation = standarddeviation(list)
    integratedautocorrelationvalues = [1/2]
    for lag in range(1, len(list)):
        normalisedautocorrelation = autocorrelation(list, lag)/selfcorrelation
        integratedautocorrelationvalues.append(integratedautocorrelationvalues[-1] + normalisedautocorrelation)

    return integratedautocorrelationvalues

def error(list):
    baseerror = standarddeviation(list)/((len(list))**0.5)
    integratedautocorrelationvalues = integratedautocorrelation(list)
    for i in range(1,len(integratedautocorrelationvalues)):
        if integratedautocorrelationvalues[i-1] > integratedautocorrelationvalues[i]:
            newerror = ((2*integratedautocorrelationvalues[i])**0.5)*baseerror
            break
    else:
        newerror = ((2*integratedautocorrelationvalues[-1])**0.5)*baseerror
    
    return newerror

def plotwitherror(xvals, yvals, err, title, filename):
    _, ax = plt.subplots()
    gfg = make_interp_spline(xvals, yvals, k=3) 
    xnew = np.linspace(min(xvals), max(xvals), 300) 
    ynew = gfg(xnew)
    ax.plot(xnew, ynew, color='blue', linestyle='-.')
    ax.errorbar(xvals, yvals, yerr=err, fmt='.', color='blue', capsize=2)
    ax.set_title(title)
    plt.savefig(filename, dpi=500)
    plt.show()

def read_16sites(file):
    with open(file, 'r') as f:
        data = f.read()
    
    rows = data.split('\n')[1:-1]
    elements = [row.split(' ')[1:-9] for row in rows]
    corrs = [[] for e in elements[1]]

    for element in elements:
        i=0
        for entry in element:
            corrs[i].append(float(entry))
            i+=1

    averages = []
    for corr in corrs:
        averages.append(average(corr))
    
    return averages

def read_32sites(file):
    with open(file, 'r') as f:
        data = f.read()
    
    rows = data.split('\n')[1:-1]
    elements = [row.split(' ')[1:-17] for row in rows]
    corrs = [[] for e in elements[1]]

    for element in elements:
        i=0
        for entry in element:
            corrs[i].append(float(entry))
            i+=1

    averages = []
    for corr in corrs:
        averages.append(average(corr))
    
    return averages