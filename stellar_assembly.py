import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import stats
import plotmedian as pm
# start up your interface of choice and define a helper function, whose purpose is to make a HTTP GET request to a specified URL ("endpoint"), and verify that the response is successful.
#---------------------------------------------------------------------------------------------
import requests

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"f0f388f8e939bf305fa9346ce99cecd1"}

def get(path, params=None):
    # make HTTP GET request to path
    headers = {"api-key":"f0f388f8e939bf305fa9346ce99cecd1"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

redshift = np.array([1.00, 0.95, 0.92, 0.89, 0.85, 0.82, 0.79, 0.76, 0.73, 0.70, 0.68, 0.64, 0.62, 0.60,
                     0.58, 0.55, 0.52, 0.50, 0.48, 0.46, 0.44, 0.42, 0.40, 0.38, 0.36, 0.35, 0.33, 0.31, 0.30,
                     0.27, 0.26, 0.24, 0.23, 0.21, 0.20, 0.18, 0.17, 0.15, 0.14, 0.13, 0.11, 0.10, 0.08, 0.07,
                     0.06, 0.05, 0.03, 0.02, 0.01, 0.00])

# cmap_reversed = matplotlib.cm.get_cmap('jet')
#2Dplot (star):   ( his2d + contour)
#---------------------------------------------------------------------------------------------

snap = 99
output = '/home/mraouf/Dropbox/KASI/TNG_Project/Raouf/plots/VYGs/'
#baseUrl = "http://www.tng-project.org/api/TNG300-1/"
stellar_file= '/media/mraouf/mraouf3/stellar_assembly.hdf5'
f = h5py.File(stellar_file, "r")                                       
print("Keys: %s" % f.keys())


StellarMassAfterInfall = f['Snapshot_99']['StellarMassAfterInfall'].value
StellarMassBeforeInfall= f['Snapshot_99']['StellarMassBeforeInfall'].value
StellarMassExSitu= f['Snapshot_99']['StellarMassExSitu'].value
StellarMassFormedOutsideGalaxies= f['Snapshot_99']['StellarMassFormedOutsideGalaxies'].value
StellarMassFromCompletedMergers= f['Snapshot_99']['StellarMassFromCompletedMergers'].value
StellarMassFromCompletedMergersMajor= f['Snapshot_99']['StellarMassFromCompletedMergersMajor'].value
StellarMassFromCompletedMergersMajorMinor= f['Snapshot_99']['StellarMassFromCompletedMergersMajorMinor'].value
StellarMassFromFlybys= f['Snapshot_99']['StellarMassFromFlybys'].value
StellarMassFromFlybysMajor= f['Snapshot_99']['StellarMassFromFlybysMajor'].value
StellarMassFromFlybysMajorMinor= f['Snapshot_99']['StellarMassFromFlybysMajorMinor'].value
StellarMassFromOngoingMergers= f['Snapshot_99']['StellarMassFromOngoingMergers'].value
StellarMassFromOngoingMergersMajor= f['Snapshot_99']['StellarMassFromOngoingMergersMajor'].value
StellarMassFromOngoingMergersMajorMinor= f['Snapshot_99']['StellarMassFromOngoingMergersMajorMinor'].value
StellarMassInSitu= f['Snapshot_99']['StellarMassInSitu'].value
StellarMassTotal= f['Snapshot_99']['StellarMassTotal'].value



f0, ((ax1)) = plt.subplots(1, 1,figsize=(11,10))


#InSitu
#*************************
ax1 =plt.subplot2grid((2, 2), (0, 0), colspan=1,rowspan=1)

w = np.where(StellarMassTotal > 0)[0]
X = np.log10(StellarMassTotal[w] *1e10)
Y = np.log10(StellarMassInSitu[w] *1e10)
xrange = [X.min(),X.max()]
yrange = [4,Y.max()]

#ax1.scatter(X,Y, marker='^', s=10, c = 'grey', alpha=0.5)
counts,ybins,xbins,image = ax1.hist2d(X,Y,bins=[100,100],cmin = 2,norm=LogNorm(),cmap = 'jet',alpha = 0.85,range=[xrange,yrange])
pm.plotmedian(X,Y,c='r',ax=ax1,bins=8,stat='median')


plt.ylabel(r'$\mathrm{M_{star}\ InSitu}$', size=15)  # Set the y...
plt.xlabel(r'$\log_{10}\ \mathrm{M_{*}\ [M_{\odot}]}$', size=15)  # and the x-axis labels

ax2 =plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)


Y = np.log10(StellarMassExSitu[w] *1e10)
yrange = [4,Y.max()]
#ax2.scatter(X,Y, marker='^', s=10, c = 'grey', alpha=0.5)
counts,ybins,xbins,image = ax2.hist2d(X,Y,bins=[100,100],cmin = 2 ,norm=LogNorm(),cmap = 'jet',alpha = 0.85,range=[xrange,yrange])

pm.plotmedian(X,Y,c='r',ax=ax2,bins=8,stat='median')



plt.ylabel(r'$\mathrm{M_{star}\ ExSitu}$', size=15)  # Set the y...
plt.xlabel(r'$\log_{10}\ \mathrm{M_{*}\ [M_{\odot}]}$', size=15)  # and the x-axis labels

ax3 =plt.subplot2grid((2, 2), (1, 0), colspan=1,rowspan=1)

#Y = np.log10(StellarMassFromFlybys[w] *1e10)
#Y = np.log10(StellarMassExSitu[w] *1e10) - X
Y = StellarMassExSitu[w] / StellarMassTotal[w]
#pm.plotmedian(X,Y,c='b',ax=ax3,bins=8,stat='median')


#Y1 = f['Snapshot_93']['StellarMassExSitu'].value / f['Snapshot_93']['StellarMassTotal'].value
#X1 = f['Snapshot_93']['StellarMassTotal'].value
#w1 = np.where(X1 > 0)[0]
#pm.plotmedian(np.log10(X1[w1]*1e10),(Y1[w1]),c='r',ax=ax3,bins=8,stat='median')

yrange = [0,1]
#ax3.scatter(X,Y, marker='^', s=10, c = 'grey', alpha=0.5)
counts,ybins,xbins,image = ax3.hist2d(X,Y,bins=[100,100],cmin = 2,norm=LogNorm(),cmap = 'jet',alpha = 0.85,range=[xrange,yrange])


plt.ylabel(r'$\mathrm{M_{star}\ Ex\ Situ\ Fraction}$', size=15)  # Set the y...
plt.xlabel(r'$\log_{10}\ \mathrm{M_{*}\ [M_{\odot}]}$', size=15)  # and the x-axis labels

ax4 =plt.subplot2grid((2, 2), (1, 1), colspan=1,rowspan=1)

X = np.log10(StellarMassInSitu[w] *1e10)
Y = np.log10(StellarMassExSitu[w] *1e10)
yrange = [5.,Y.max()]
xrange = [5.,Y.max()]
#ax4.scatter(X,Y, marker='^', s=10, c = 'grey', alpha=0.5)
counts,ybins,xbins,image = ax4.hist2d(X,Y,bins=[100,100],cmin = 2,norm=LogNorm(),cmap = 'jet',alpha = 0.85,range=[xrange,yrange])

#xx = [2,2,2]
#yy = [0,1,3]
#plt.plot(xx,yy,'k:', lw = 2)
plt.plot(xrange,yrange,'k--', lw = 2)


#pm.plotmedian(X,Y,c='r',ax=ax4,bins=8,stat='median')


plt.ylabel(r'$\mathrm{Stellar\ Mass\ In\ Situ}$', size=15)  # Set the y...
plt.xlabel(r'$\log_{10}\ \mathrm{Stellar\ Mass\ Ex\ Situ}$', size=15)  # and the x-axis labels


#leg = plt.legend(loc='upper left',labelspacing = 0.05)
#leg.draw_frame(False)  # Don't want a box frame
#for t in leg.get_texts():  # Reduce the size of the text
    #t.set_fontsize('x-large')

#ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)

filename = output+'VYGs_Mstar_InSitu'
plt.savefig(filename + '.pdf')
plt.close()

