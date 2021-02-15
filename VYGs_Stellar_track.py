import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import stats
import plotmedian as pm
import illustris_python as il
from matplotlib import colors, ticker, cm
import random
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

red_z = np.array([1.00, 0.95, 0.92, 0.89, 0.85, 0.82, 0.79, 0.76, 0.73, 0.70, 0.68, 0.64, 0.62, 0.60,
                     0.58, 0.55, 0.52, 0.50, 0.48, 0.46, 0.44, 0.42, 0.40, 0.38, 0.36, 0.35, 0.33, 0.31, 0.30,
                     0.27, 0.26, 0.24, 0.23, 0.21, 0.20, 0.18, 0.17, 0.15, 0.14, 0.13, 0.11, 0.10, 0.08, 0.07,
                     0.06, 0.05, 0.03, 0.02, 0.01, 0.00])


# Read data
# ---------------------------------------------------------------------------
output = '/home/mraouf/Dropbox/KASI/TNG_Project/Raouf/plots/VYGs/'
#data_path ='/home/mraouf/Dropbox/KASI/TNG_Project/Salimi/data/TNG_50/group_properties_TNG100_2.csv' 
#basePath  ='/home/mraouf/KASI_Projects/work_space/Illustris/illustris_data/ILLUSTRIS_TNG/TNG100-1/output/'
basePath  ='/home/mraouf/KASI_Projects/work_space/Illustris/illustris_data/ILLUSTRIS_TNG/TNG50-1/output/'
#output_data = '/home/mraouf/Dropbox/KASI/TNG_Project/Raouf/galaxies_properties_all_TNG100.asc'
output_data = '/home/mraouf/Dropbox/KASI/TNG_Project/Raouf/galaxies_properties_all_TNG50.asc'

#group_properties_TNG100_VYGs.txt
# write your data_path in your system
#groupdata_TNG = pd.read_csv(data_path)


#general conditions:
# ---------------------------------------------------------------------------

Group = il.groupcat.loadHalos(basePath,99,fields=['GroupFirstSub','Group_M_Crit200'])
Subhalo = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloGrNr','SubhaloMassInHalfRad'])
GroupFirstSub = Group['GroupFirstSub'][(Group['Group_M_Crit200']>0.01)]
#ids = random.sample(range(20000), 2000)

mass = Subhalo['SubhaloMassInHalfRad']
w = np.where(np.log10(mass*1e10) > 8)[0]

start = 0
end = len(Subhalo)
# get information of ID from URL
#for i, id in enumerate(range(start,start+5)) :
color = iter(cm.rainbow(np.linspace(0, 1, 5)))
#ratio = 1.0/5.0
rassembly = []
rsf = []
j = 0
halfmassrad= [];halfmassrad_stars= [];len_stars= [];mass_stars= [];sfr= [];sfrinhalfrad= [];sfrinmaxrad= [];sfrinrad= [];stellarphotometrics_u= [];stellarphotometrics_b= [];stellarphotometrics_v= [];stellarphotometrics_k= [];stellarphotometrics_g= [];stellarphotometrics_r= [];stellarphotometrics_i= [];stellarphotometrics_z= [];massinhalfrad_stars=[];id1= []
jj = 0
#for i in GroupFirstSub:
#for i in range(start,end) :
for i in w:
 #if i > 32255:   
  #id  = GroupFirstSub[i]
  fields1 = ['SubhaloID','NextProgenitorID','MainLeafProgenitorID','FirstProgenitorID','SubhaloMassType','SnapNum','SubfindID','Group_M_Mean200','SubhaloGrNr','SubhaloSFRinRad','SubhaloHalfmassRadType','SubhaloMassInRad','SubhaloMassInRadType','SubhaloMassInHalfRad']
  tree = il.sublink.loadTree(basePath,99,i,fields=fields1, onlyMPB=True)
  tree_all = il.sublink.loadTree(basePath,99,i,fields=fields1, onlyMPB=False)
  #ptNumStars = il.snapshot.partTypeNum('stars') # 4
  #stars_mass = tree['SubhaloMassType'][ptNumStars]
  if tree != None:
    mass_stellar = tree['SubhaloMassInRad']
    massinhalfrad = tree['SubhaloMassInHalfRad']
    stars_mass = tree['SubhaloMassInRadType'] 
    stars_sf = tree['SubhaloSFRinRad']
    flag_93 = np.where(tree['SnapNum'] == 93)
    flag_93_99 = np.where(tree['SnapNum'] >= 93)
    flag_99 = np.where(tree['SnapNum'] == 99)
    
    stars_mass_all = tree_all['SubhaloMassInRadType'] 
    #stars_sf_all = tree_all['SubhaloSFRinRad']
    #flag_93_all = np.where(tree['SnapNum'] == 93)
    
    #flag_94_all = np.where(tree_all['SnapNum'] == 94)
    #flag_95_all = np.where(tree_all['SnapNum'] == 95)
    #flag_96_all = np.where(tree_all['SnapNum'] == 96)
    #flag_97_all = np.where(tree_all['SnapNum'] == 97)
    #flag_98_all = np.where(tree_all['SnapNum'] == 98)
    #flag_99_all = np.where(tree_all['SnapNum'] == 99)
    
   #if len() 
    sf_93  = stars_sf[flag_93_99]
    lookback_time = [0.0,0.136,0.340,0.475,0.676,0.810,1.008]
    #mass_sf_93 = sum(stars_mass_all[flag_93_all])
    mass_93  = stars_mass[flag_93]
    mass_99 = stars_mass[flag_99]
    MS = np.log10(mass_stellar[flag_99]* 1e10)
    MS_massinhalfrad = np.log10(massinhalfrad[flag_99]* 1e10)
    if ((len(flag_93[0]) > 0)):
        ki = (mass_93[0][4] *1e10) + sum([stars_sf[j] * ((lookback_time[j+1]-lookback_time[j])* 1e9)  for j in range(0,len(stars_sf[flag_93_99])-1)])

        sf_93  = stars_sf[flag_93]
        sf_99 = stars_sf[flag_99]
    
        r_assembly_dir = mass_99 / mass_93
        #r_sf_dir  = mass_99 / mass_sf_93
        r_sf = (mass_99 * 1e10) / ki
        rassembly.append(r_assembly_dir[0][4])
        id1.append(i)
        rsf.append(r_sf[0][4])
        mass_stars.append(MS)
        massinhalfrad_stars.append(MS_massinhalfrad)

        jj = jj+1
        print jj, i, r_assembly_dir[0][4], r_sf[0][4], MS_massinhalfrad
      
cc =  np.array([id1,rassembly,rsf,mass_stars,massinhalfrad_stars]).T
np.savetxt(output_data, cc, fmt='%i %1.8f %1.8f  %1.8f %1.8f', delimiter='\t',header='id\trassembly\trsf\tmass_stars\tmassinhalfrad_stars')

r_assembly = np.array(rassembly)
r_sf       = np.array(rsf)
xrange = [0,3]
yrange = [0,3]
plt.scatter(r_assembly,r_sf,c = 'grey', s = 3)

xx = [2,2,2]
yy = [0,1,3]
plt.plot(xx,yy,'k:', lw = 2)
plt.plot(xrange,yrange,'k--', lw = 2)

plt .axis([xrange[0],xrange[1], yrange[0],yrange[1]])
plt.xlabel(r'$\mathrm{r_{assembly}}$', size=15)  
plt.ylabel(r'$\mathrm{r_{SF}}$', size=15)  # and the y-axis labels

#plt.plot(redshift,star_mass_history,c = c,label="z = 0")


#plt.xlabel(r'$\mathrm{z}$', size=15)  
#plt.ylabel(r'$\log_{10}\ \mathrm{M_{*}\ [M_{\odot}]}$', size=15)  # and the y-axis labels


#leg = plt.legend(loc='upper left',labelspacing = 0.05)
#leg.draw_frame(False)  # Don't want a box frame
#for t in leg.get_texts():  # Reduce the size of the text
    #t.set_fontsize('x-large')

#ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1,rowspan=1)

filename = output+'VYGs_rass_rsf_z_all_TNG50'
plt.savefig(filename + '.pdf')
plt.close()

'''

  url = 'https://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/' + str(id)
  sub = get(url) # get json response of subhalo properties
  if (sub['stellarphotometrics_i'] <= -13):
    j = j + 1  
    print j, i, id
    # prepare dict to hold result arrays
    fields = ['snap','id','mass_gas','mass_stars','mass_dm','mass_bhs','sfr','spin_x','spin_y','spin_z',
              'starmetallicity','gasmetallicity','mass_log_msun','massinhalfrad_stars','stellarphotometrics_i','sfrinrad']
    r = {}
    for field in fields:
        r[field] = []
    ii = 0
    redshift = []
    #c = next(color)
    while sub['prog_snap'] >= 92:
        redshift.append(red_z[49 - ii])
        print 'start', sub['desc_snap'],red_z[49 - ii]
        ii = ii + 1
        for field in fields:
            r[field].append(sub[field])
        # request the full subhalo details of the descendant by following the sublink URL
        sub = get(sub['related']['sublink_progenitor'])

    gas_mass_history  =  (np.array(r['mass_gas']) * 1e10/0.704)
    star_mass_history =  (np.array(r['mass_stars']) * 1e10/0.704)
    sfr_history = (np.array(r['sfrinrad']))
    #redshift = np.array(redshift)
    if (len(star_mass_history) > 6):
      rassembly1 = star_mass_history[0]/ star_mass_history[6]
      
      rsf1 = sum((star_mass_history/star_mass_history[6]) -1)
      rsfr = sfr_history[0]/sfr_history[6]
      print 'r_main, rsf, rsfr',rassembly1, rsf1, r_assembly_dir[0][4]
      rassembly.append(rassembly1)
      id1.append(id)
      rsf.append(rsf1)     
      halfmassrad.append( sub['halfmassrad']  );halfmassrad_stars.append( sub['halfmassrad_stars']  );len_stars.append( sub['len_stars']  );mass_stars.append( sub['mass_stars']  );sfr.append( sub['sfr']  );sfrinhalfrad.append( sub['sfrinhalfrad']  );sfrinmaxrad.append(  sub['sfrinmaxrad'] );sfrinrad.append(sub['sfrinrad']   );stellarphotometrics_u.append( sub['stellarphotometrics_u']  );stellarphotometrics_b.append( sub['stellarphotometrics_b']  );stellarphotometrics_v.append(  sub['stellarphotometrics_v'] );stellarphotometrics_k.append(  sub['stellarphotometrics_k'] );stellarphotometrics_g.append(  sub['stellarphotometrics_g'] );stellarphotometrics_r.append(sub['stellarphotometrics_r']   );stellarphotometrics_i.append( sub['stellarphotometrics_i']  );stellarphotometrics_z.append(  sub['stellarphotometrics_z'] );massinhalfrad_stars.append(sub['massinhalfrad_stars'])

    
cc =  np.array([id1,rassembly,rsf,halfmassrad,halfmassrad_stars,len_stars,mass_stars,sfr,sfrinhalfrad,sfrinmaxrad,sfrinrad,stellarphotometrics_u,stellarphotometrics_b,stellarphotometrics_v,stellarphotometrics_k,stellarphotometrics_g,stellarphotometrics_r,stellarphotometrics_i,stellarphotometrics_z,massinhalfrad_stars]).T
np.savetxt(output_data, cc, fmt='%i %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f', delimiter='\t',header='id\trassembly\trsf\thalfmassrad\thalfmassrad_stars\tlen_stars\tmass_stars\tsfr\tsfrinhalfrad\tsfrinmaxrad\tsfrinrad\tstellarphotometrics_u\tstellarphotometrics_b\tstellarphotometrics_v\tstellarphotometrics_k\tstellarphotometrics_g\tstellarphotometrics_r\tstellarphotometrics_i\tstellarphotometrics_z\tmassinhalfrad_stars')
'''
