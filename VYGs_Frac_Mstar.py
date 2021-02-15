# -*- coding: utf-8 -*-
import matplotlib
# import pylab as plt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.ndimage
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
# from matplotlib.mlab import bivariate_normal
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter
from astropy.io import fits
import plotmedian as pm
import pandas as pd
from io import StringIO
from pandas.plotting import scatter_matrix
import astropy.stats as ass
import h5py
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
def customizedLogFormat(x,pos):
  decimalplaces = int(np.maximum(-np.log10(x),0))
  formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
  return formatstring.format(x)



from scipy.stats import beta
def binom_interval(success, total, confint=0.95):
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total - success + 1)
    upper = beta.ppf(1 - quantile, success + 1, total - success)
    return (lower, upper)

def read_data(filename):

    data = np.genfromtxt(filename,names=True,dtype=None) #dtype=None to be able to read mixed data types form the file
    return data


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


VYGs_data_cent = '/home/mraouf/Dropbox/KASI/VYGs_Project/galaxies_properties_cent_TNG100.asc'

VYGs_data_sat = '/home/mraouf/Dropbox/KASI/VYGs_Project/galaxies_properties_sat_TNG100.asc'
VYGs_data_all = '/home/mraouf/Dropbox/KASI/VYGs_Project/galaxies_properties_all_TNG100.asc'
VYGs_TNG50_all = '/home/mraouf/Dropbox/KASI/VYGs_Project/galaxies_properties_all_TNG50.asc'


#id rassembly rsf mass_stars massinhalfrad_stars
#df = pd.read_csv(VYGs_data, header=0)
df = read_data(VYGs_data_cent)
df_sat = read_data(VYGs_data_sat)
df_all = read_data(VYGs_data_all)
df_TNG50 =read_data(VYGs_TNG50_all)
output = './plots/VYGs/'

# *************************************
#Fraction of Mis-align vs BGG mass_bin in different mass
# *************************************
f, ((ax1)) = plt.subplots(1, 1,figsize=(6,5))
bin_ms = 5

Q = np.where((df['rsf'] <= 100.)&(df['mass_stars'] > 8.5))
Q_sat = np.where((df_sat['rsf'] <= 100.)&(df_sat['mass_stars'] > 8.5))
Q_all = np.where((df_all['rsf'] <= 100.)&(df_all['mass_stars'] > 8.5))


Q_TNG50 = np.where((df_TNG50['rsf'] <= 100.)&(df_TNG50['mass_stars'] > 8.5))
#ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=1,rowspan=1)


VYG = df['rsf'][Q]
Mass = (df['mass_stars'][Q])

VYG_sat = df_sat['rsf'][Q_sat]
Mass_sat = (df_sat['mass_stars'][Q_sat])

VYG_all = df_all['rsf'][Q_all]
Mass_all= (df_all['mass_stars'][Q_all])

VYG_TNG50 = df_TNG50['rsf'][Q_TNG50]
Mass_TNG50 = (df_TNG50['mass_stars'][Q_TNG50])

# dist12 =
MinRange = 8.5
MaxRange = 11.
Interval = 0.5
Nbins = int((MaxRange-MinRange)/Interval)
Range = np.arange(MinRange, MaxRange, Interval)

xrange = [8.5,10.5]
yrange = [0.0001,0.025]
#color = iter(cm.rainbow(np.linspace(0, 1, 5)))
#c = next(color)
#*******************************************
mass_bin = []
fraction_VYGs = []
fraction_NonVY = []
bino_error = []

fraction_VYGs_sat = []
fraction_NonVY_sat = []
bino_error_sat = []

fraction_VYGs_all = []
fraction_NonVY_all = []
bino_error_all = []


fraction_VYGs_TNG50 = []
fraction_NonVY_TNG50 = []
bino_error_TNG50 = []



for i in range(Nbins-1):

    w = np.where((Mass >= Range[i]) & (Mass < Range[i+1]))[0]
    w_sat = np.where((Mass_sat >= Range[i]) & (Mass_sat < Range[i+1]))[0]
    w_all = np.where((Mass_all >= Range[i]) & (Mass_all < Range[i+1]))[0]
    w_TNG50 = np.where((Mass_TNG50 >= Range[i]) & (Mass_TNG50 < Range[i+1]))[0]

    if len(w) > 0:
        wNonVY = np.where((Mass >= Range[i]) & (Mass < Range[i+1]) & (VYG < 2))[0]
        wVYG = np.where((Mass >= Range[i]) & (Mass < Range[i+1]) & (VYG > 2.))[0]
        Frac_NonVY = 1.0*len(wNonVY) / np.float(len(w))
        Frac_VYG = 1.0*len(wVYG) / np.float(len(w))
        #print len(wNonVY), len(wVYG), Frac_NonVY, Frac_VYG, len(w)
        fraction_VYGs.append(Frac_VYG)
        fraction_NonVY.append(Frac_NonVY)
        # bino_error.append(ass.binom_conf_interval(Frac_VYG, len(w), 0.68269, interval='jeffreys'))
        bino_error.append(np.sqrt(Frac_VYG*(1-Frac_VYG)/len(w)))
    else:
        bino_error.append(0.0)
        fraction_NonVY.append(0.0)
        fraction_VYGs.append(0.0)

    if len(w_sat) > 0:
        wNonVY_sat = np.where((Mass_sat >= Range[i]) & (Mass_sat < Range[i+1]) & (VYG_sat < 2))[0]
        wVYG_sat = np.where((Mass_sat >= Range[i]) & (Mass_sat < Range[i+1]) & (VYG_sat > 2.))[0]
        Frac_NonVY_sat = 1.0*len(wNonVY_sat) / np.float(len(w_sat))
        Frac_VYG_sat= 1.0*len(wVYG_sat) / np.float(len(w_sat))
        #print len(wNonVY), len(wVYG), Frac_NonVY, Frac_VYG, len(w)
        fraction_VYGs_sat.append(Frac_VYG_sat)
        fraction_NonVY_sat.append(Frac_NonVY_sat)
        # bino_error_sat.append(ass.binom_conf_interval(Frac_VYG_sat, len(w), 0.68269, interval='jeffreys'))
        bino_error_sat.append(np.sqrt(Frac_VYG_sat*(1-Frac_VYG_sat)/len(w_sat)))
    else:
        bino_error_sat.append(0.0)
        fraction_NonVY_sat.append(0.0)
        fraction_VYGs_sat.append(0.0)
    if len(w_all) > 0:
        wNonVY_all = np.where((Mass_all >= Range[i]) & (Mass_all < Range[i+1]) & (VYG_all < 2))[0]
        wVYG_all = np.where((Mass_all >= Range[i]) & (Mass_all < Range[i+1]) & (VYG_all > 2.))[0]
        Frac_NonVY_all = 1.0*len(wNonVY_all) / np.float(len(w_all))
        Frac_VYG_all= 1.0*len(wVYG_all) / np.float(len(w_all))
        #print len(wNonVY), len(wVYG), Frac_NonVY, Frac_VYG, len(w)
        fraction_VYGs_all.append(Frac_VYG_all)
        fraction_NonVY_all.append(Frac_NonVY_all)
        # bino_error_all.append(ass.binom_conf_interval(Frac_VYG_all, len(w), 0.68269, interval='jeffreys'))
        bino_error_all.append(np.sqrt(Frac_VYG_all*(1-Frac_VYG_all)/len(w_all)))
    else:
        bino_error_all.append(0.0)
        fraction_NonVY_all.append(0.0)
        fraction_VYGs_all.append(0.0)

    if len(w_TNG50) > 0:
        wNonVY_TNG50 = np.where((Mass_TNG50 >= Range[i]) & (Mass_TNG50 < Range[i+1]) & (VYG_TNG50 < 2))[0]
        wVYG_TNG50 = np.where((Mass_TNG50 >= Range[i]) & (Mass_TNG50 < Range[i+1]) & (VYG_TNG50 > 2.))[0]
        Frac_NonVY_TNG50 = 1.0*len(wNonVY_TNG50) / np.float(len(w_TNG50))
        Frac_VYG_TNG50 = 1.0*len(wVYG_TNG50) / np.float(len(w_TNG50))
        #print len(wNonVY), len(wVYG), Frac_NonVY, Frac_VYG, len(w)
        fraction_VYGs_TNG50.append(Frac_VYG_TNG50)
        fraction_NonVY_TNG50.append(Frac_NonVY_TNG50)
        # bino_error_TNG50.append(ass.binom_conf_interval(Frac_VYG_TNG50, len(w), 0.68269, interval='jeffreys'))
        bino_error_TNG50.append(np.sqrt(Frac_VYG_TNG50*(1-Frac_VYG_TNG50)/len(w_TNG50)))
    else:
        bino_error_TNG50.append(0.0)
        fraction_NonVY_TNG50.append(0.0)
        fraction_VYGs_TNG50.append(0.0)
    mass_bin.append((Range[i] + Range[i+1]) / 2.0)
            # print '  ', mass_bin[i], Fraction[i]


mass_bin = np.array(mass_bin)
fraction_VYGs = np.array(fraction_VYGs)
fraction_NonVY = np.array(fraction_NonVY)
yy_err = np.array(bino_error)

fraction_VYGs_sat = np.array(fraction_VYGs_sat)
fraction_NonVY_sat = np.array(fraction_NonVY_sat)
yy_err_sat = np.array(bino_error_sat)

fraction_VYGs_all = np.array(fraction_VYGs_all)
fraction_NonVY_all = np.array(fraction_NonVY_all)
yy_err_all = np.array(bino_error_all)

mass_bin = np.array(mass_bin)
fraction_VYGs_TNG50 = np.array(fraction_VYGs_TNG50)
fraction_NonVY_TNG50 = np.array(fraction_NonVY_TNG50)
yy_err_TNG50 = np.array(bino_error_TNG50)

w = np.where(fraction_VYGs > 0)[0]
plt.plot(mass_bin[w], fraction_VYGs[w], color = 'blue', label = 'TNG-100(Central)')
yvalU = fraction_VYGs[w]+(yy_err[w])
yvalL = fraction_VYGs[w]-(yy_err[w])
plt.fill_between(mass_bin[w], yvalU,yvalL, facecolor='blue', alpha=0.15)


w = np.where(fraction_VYGs_sat > 0)[0]
plt.plot(mass_bin[w], fraction_VYGs_sat[w], color = 'red', label = 'TNG-100(Satelite)')
yvalU = fraction_VYGs_sat[w]+(yy_err_sat[w])
yvalL = fraction_VYGs_sat[w]-(yy_err_sat[w])
plt.fill_between(mass_bin[w], yvalU,yvalL, facecolor='red', alpha=0.15)

w = np.where(fraction_VYGs_all > 0)[0]
plt.plot(mass_bin[w], fraction_VYGs_all[w], color = 'black', label = 'TNG-100(all)')
yvalU = fraction_VYGs_all[w]+(yy_err_all[w])
yvalL = fraction_VYGs_all[w]-(yy_err_all[w])
plt.fill_between(mass_bin[w], yvalU,yvalL, facecolor='grey', alpha=0.15)


w = np.where(fraction_VYGs_TNG50 > 0)[0]
plt.plot(mass_bin[w], fraction_VYGs_TNG50[w], color = 'green', label = 'TNG-50(All)')
yvalU = fraction_VYGs_TNG50[w]+(yy_err_TNG50[w])
yvalL = fraction_VYGs_TNG50[w]-(yy_err_TNG50[w])
plt.fill_between(mass_bin[w], yvalU,yvalL, facecolor='green', alpha=0.15)


#plt.text(8.2,0.023, 'TNG-100(All Galaxies)', fontsize=15, rotation=0)

plt.yscale("log")
# plt.xscale("log")
plt.axis([xrange[0],xrange[1], yrange[0],yrange[1]])

plt.ylabel(r'$\mathrm{Fraction\ of\ VYGs}$', size=15)  # Set the y...
plt.xlabel(r'$\mathrm{M_{star}\ [M_{\odot}]}$', size=15)  # and the x-axis labels
leg = plt.legend(loc='top right',labelspacing = 0.05)
leg.draw_frame(False)  # Don't want a box frame
for t in leg.get_texts():  # Reduce the size of the text
    t.set_fontsize('large')

plt.subplots_adjust(wspace=0.1,hspace=0.35, bottom=0.12, right=0.98,left=0.15, top=0.98)


# Save to a File
filename = output+'VYGs-Frac-Mass'
plt.savefig(filename + '.pdf')
plt.close()

# *******************************************************************************************
f, ((ax1)) = plt.subplots(1, 1,figsize=(6,5))
X  = df_all['rsf']
Y  = df_all['rassembly']
Z  =df_all['mass_stars']
xrange = [0,7]
yrange = [0,7]
# plt.scatter(r_assembly,r_sf,c = 'grey', s = 3)
# counts,ybins,xbins,image = plt.hist2d(X,Y ,bins=[70,70],norm=LogNorm(),cmin = 1,cmap = 'jet',alpha = 0.85, range=[xrange,yrange])
plt.hexbin(X,Y, C=Z, gridsize=200, cmap='jet', bins=None,edgecolors=None,mincnt = 0.1)
cb = plt.colorbar()
cb.set_label(r'$\log_{10}\ (\mathrm{M_{star}\ [M_{\odot}]})$', rotation=270,labelpad=20)
xx = [2,2,2]
yy = [0,1,8]
plt.plot(xx,yy,'k:', lw = 2)
plt.plot(xrange,yrange,'k--', lw = 2)

plt .axis([xrange[0],xrange[1], yrange[0],yrange[1]])
plt.ylabel(r'$\mathrm{Assembly\ growth}}$', size=15)
plt.xlabel(r'$\mathrm{Star\ formation\ growth}$', size=15)  # and the y-axis labels

#plt.plot(redshift,star_mass_history,c = c,label="z = 0")


#plt.xlabel(r'$\mathrm{z}$', size=15)
#plt.ylabel(r'$\log_{10}\ \mathrm{M_{*}\ [M_{\odot}]}$', size=15)  # and the y-axis labels


#leg = plt.legend(loc='upper left',labelspacing = 0.05)
#leg.draw_frame(False)  # Don't want a box frame
#for t in leg.get_texts():  # Reduce the size of the text
    #t.set_fontsize('x-large')

#ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)

filename = output+'VYGs_rass_rsf_TNG100_all'
plt.savefig(filename + '.pdf')
plt.close()

# *******************************************************************************************
f, ((ax1)) = plt.subplots(1, 1,figsize=(6,5))
X  = df_all['mass_stars']
Y  = df_all['rassembly']/df_all['rsf']
xrange = [8,12]
yrange = [1,10]
# plt.scatter(r_assembly,r_sf,c = 'grey', s = 3)
counts,ybins,xbins,image = plt.hist2d(X,Y ,bins=[70,70],norm=LogNorm(),cmin = 1,cmap = 'jet',alpha = 0.85, range=[xrange,yrange])
# plt.hexbin(X,Y, C=dvx_star, gridsize=100, cmap='jet_r', bins=None,edgecolors=None,mincnt = 0.1)
# cb = plt.colorbar()

plt .axis([xrange[0],xrange[1], yrange[0],yrange[1]])
plt.ylabel(r'$\mathrm{Assembly/SF\  growth}}$', size=15)
plt.xlabel(r'$\mathrm{M_{star} [M_{\odot}]}$', size=15)  # and the y-axis labels

#plt.plot(redshift,star_mass_history,c = c,label="z = 0")


#plt.xlabel(r'$\mathrm{z}$', size=15)
#plt.ylabel(r'$\log_{10}\ \mathrm{M_{*}\ [M_{\odot}]}$', size=15)  # and the y-axis labels


#leg = plt.legend(loc='upper left',labelspacing = 0.05)
#leg.draw_frame(False)  # Don't want a box frame
#for t in leg.get_texts():  # Reduce the size of the text
    #t.set_fontsize('x-large')

#ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)

filename = output+'VYGs_rass_rsf_mass_TNG100_all'
plt.savefig(filename + '.pdf')
plt.close()





'''
#****************************************************************************8
f, ((ax1)) = plt.subplots(1, 1,figsize=(6,5))


id =  df['id']
MassInSitu = np.log10(StellarMassInSitu[id] *1e10)
MassExSitu = np.log10(StellarMassExSitu[id] *1e10)
r_assembly = df['rassembly']
r_sf       = df['rsf']
xrange = [0.8,6]
yrange = [0.12,0.4]
plt.scatter(r_assembly,r_sf,c = MassInSitu, s = 3, cmap = 'jet')
cb = plt.colorbar()
cb.set_label(r'$\log_{10} (\mathrm{Mass\ InSitu})$', rotation=270,labelpad=20)
xx = [2,2,2]
yy = [0,1,3]
plt.plot(xx,yy,'k:', lw = 2)
#plt.plot(xrange,yrange,'k--', lw = 2)

plt .axis([xrange[0],xrange[1], yrange[0],yrange[1]])
plt.xlabel(r'$\mathrm{r_{assembly}}$', size=15)
plt.ylabel(r'$\mathrm{r_{SF}}$', size=15)  # and the y-axis labels

#plt.text(-1.95, 0.105, 'GAMA'+str(data['CATAID'][ID_reg][0]), fontsize=8, rotation=0)


#leg = plt.legend(loc='upper left',labelspacing = 0.05)
#leg.draw_frame(False)  # Don't want a box frame
#for t in leg.get_texts():  # Reduce the size of the text
    #t.set_fontsize('x-large')

#ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)

filename = output+'VYGs_rass_rsf_InSitu'
plt.savefig(filename + '.pdf')
plt.close()


#****************************************************************************8
f, ((ax1)) = plt.subplots(1, 1,figsize=(6,5))


id =  df['id']
#MassInSitu = np.log10(StellarMassInSitu[id] *1e10)
MassExSitu = np.log10(StellarMassExSitu[id] *1e10)
r_assembly = df['rassembly']
r_sf       = df['rsf']
xrange = [0.8,6]
yrange = [0.12,0.4]
plt.scatter(r_assembly,r_sf,c = MassExSitu, s = 3, cmap = 'jet')
cb = plt.colorbar()
cb.set_label(r'$\log_{10} (\mathrm{Mass\ ExSitu})$', rotation=270,labelpad=20)
xx = [2,2,2]
yy = [0,1,3]
plt.plot(xx,yy,'k:', lw = 2)
#plt.plot(xrange,yrange,'k--', lw = 2)

plt .axis([xrange[0],xrange[1], yrange[0],yrange[1]])
plt.xlabel(r'$\mathrm{r_{assembly}}$', size=15)
plt.ylabel(r'$\mathrm{r_{SF}}$', size=15)  # and the y-axis labels

#plt.text(-1.95, 0.105, 'GAMA'+str(data['CATAID'][ID_reg][0]), fontsize=8, rotation=0)


#leg = plt.legend(loc='upper left',labelspacing = 0.05)
#leg.draw_frame(False)  # Don't want a box frame
#for t in leg.get_texts():  # Reduce the size of the text
    #t.set_fontsize('x-large')

#ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)

filename = output+'VYGs_rass_rsf_ExSitu'
plt.savefig(filename + '.pdf')
plt.close()
'''
