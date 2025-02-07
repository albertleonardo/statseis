

# which didtance function to call?


import etannd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime as UT

fname = 'hauksson.csv'
#fname='2020_2023_ncedc.csv'

df = pd.read_csv(fname,sep=' ')
print(df)
#df = df[df['Magnitude']>=4]
print(df)

# add a small random number to avoid problems with quakes at the exact same coordinates
df['Longitude'] = df['Longitude']+np.random.random(len(df['Longitude']))/10000000000
df['Latitude'] = df['Latitude']+np.random.random(len(df['Latitude']))/10000000000


lons=df['Longitude'].to_numpy()
lats=df['Latitude'].to_numpy()
deps=df['Depth'].to_numpy()
mags=df['Magnitude'].to_numpy()
times=df['DateTime'].to_list()
ids=df['EventID'].to_numpy()

"""
ref_index = 6000

ref_lat = lats[ref_index] 
ref_lon = lons[ref_index]
ref_time= times[ref_index]
ref_dep = deps[ref_index]
ref_mag = mags[ref_index]

#print(ref_lat,ref_lon,ref_time,ref_dep)
#print(UT(ref_time))

rescaled_time,rescaled_distance,eta = etannd.estimate_eta(ref_time,ref_lat,ref_lon,ref_dep,ref_mag,times,lats,lons,deps,mags)

#print(min(eta))
"""

nndistances = []
nn_T = []
nn_R = []
nn_id = []
#for i in range(len(lons)):

nndistances.append(1)
nn_T.append(1)
nn_R.append(1)
nn_id.append(1)


for i in range(1,len(lons)):
	print(i)
	ref_index = i
	ref_lat = lats[ref_index]
	ref_lon = lons[ref_index]
	ref_time= times[ref_index]
	ref_dep = deps[ref_index]
	ref_mag = mags[ref_index]

	rescaled_time,rescaled_distance,eta,t_ids = etannd.estimate_eta(ref_time,ref_lat,ref_lon,0,ref_mag,times,lats,lons,np.zeros(len(lons)),mags,ids,q=0.5)
	#nn_rescaled_time     = np.min(rescaled_time) # make it years
	#nn_rescaled_distance = np.min(rescaled_distance)
	#print(eta)
	nn_eta               = np.min(eta)
	min_eta_arg = np.argmin(eta)
	nn_rescaled_time     = rescaled_time[min_eta_arg]
	nn_rescaled_distance = rescaled_distance[min_eta_arg]
	#print(nn_eta)
	nndistances.append(nn_eta)
	#these migh not correspond to the min eta
	nn_T.append(nn_rescaled_time)
	nn_R.append(nn_rescaled_distance)
	nn_id.append(t_ids[min_eta_arg])

# add columns  to the datafram with the results and save it
# add a row to make compatible
#nndistances.append(1)
#nn_T.append(1)
#nn_R.append(1)
#nn_id.append(1)

df['nn_T'] = nn_T
df['nn_R'] = nn_R
df['nn_eta'] = nndistances
df['nn_id'] = nn_id

outname = fname.split('.')[0]+'_nnd.csv'

df.to_csv(outname)



plt.hist(np.log10(nndistances),bins=40)
plt.show()

plt.figure()
plt.hexbin(np.log10(nn_T),np.log10(nn_R), gridsize=60, cmap='jet')
plt.show()

plt.figure()
plt.scatter(np.log10(nn_T),np.log10(nn_R),alpha=0.03)
plt.show()


