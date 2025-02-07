
# Earthquake clustering using Nearest Neighbor Distance in eta space
# Copycat of Zaliapin and Ben-Zion 2013

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mbc

from obspy.core import UTCDateTime as UT

def estimate_eta(ref_time,ref_lat,ref_lon,ref_dep,ref_mag,times,lats,lons,deps,mags,ids,b=1,df=1.6,q=0.5):
	"""
	estimate rescaled times and rescaled distances, and eta=TR
	given an earthquake catalog and a single event
	"""
	#rescaled time, delays come out in days from mbc.py
	delays            = calculate_delays(times,ref_time)/365.25  # make it years
	#print(delays)
	#delays = -1*delays
	positive_delays   = delays[delays>0]
	#print(positive_delays[:20]);print(type(positive_delays))
	distances         = mbc.calculate_distances_old(lats,lons,deps,ref_lat,ref_lon,ref_dep)# ignore depths
	print(distances[:10])
	t_mags            = mags[delays>0]
	t_ids             = ids[delays>0]
	#print(t_mags)
	t_distances       = distances[delays>0]
	rescaled_time     = np.asarray(positive_delays)*10**(-q*b*t_mags)
	rescaled_distance = (t_distances**(df))*10**(-(1-q)*b*t_mags)
	eta               = rescaled_time*rescaled_distance
	# need to pass the index of the nearest neighbor for subsequent treatment
	return rescaled_time,rescaled_distance,eta,t_ids

# no worries about indexing so far
def find_nearest_neighbor(eta):
	return np.min(eta)


# replace the calcualte_delays function




def calculate_delays(times,ref_time):
	"""
	function calculate_delays
	take the inputs as strings, convert to UTCDateTime and calculate delays
	"""
	ut_times = [UT(time) for time in times]
	ref_ut_time = UT(ref_time)	

	delays = [ref_ut_time - ut_time for ut_time in ut_times]
	# delays are in seconds, convert to days
	delays = np.asarray(delays)
	delays = delays / 86400

	return delays

