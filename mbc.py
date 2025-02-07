
#############################################
# magnitud based clustering MBC             #
#############################################

import numpy as np
from obspy.core import UTCDateTime as UT

# to fix:
# use delays that imply earthquakes are only after, not before
# vectorize calculate delays


def mbc(lons,lats,deps,mags,times,ids,min_mag=3):
	"""
	function mbc: magnitude based clustering

	TO DO:
	1. check array sizes	


	:type lons: numpy.array
	:param lons:
	     longitudes of hypocenters in decimal degrees
	:type lats: numpy.array
	:param lats:
	    latitudes of hypocenters in decimal degrees
	:type deps: numpy.array
	:param deps:
	    depth of hypocenters in kilometers
	:type mags: numpy.array	
	:param mags:
	    magnitudes of earthquakes in the catalog
	    For now, it doesn't consider different scales
	:type times: list
	:param times:
	    origin times of hypocenters
	    Each entry in the list should be a string	
	:type ids: numpy.array
	:param ids:
	    Hypocenter ID's, this is optional as this 
	    routine assigns ids, if origin_ID's are available,
	    those can be used.
	:type min_mag: float
	:param min_mag:
	    Minimum magnitude to consider as mainshock
	
	.. Note::
	    The distance function is not great
	"""
	assignments = np.zeros(len(lons))
	print('initial assignments')	
	i = 1
	t_max_mag=10
	while t_max_mag > min_mag:
		assignments = one_mbc(lons,lats,deps,mags,times,ids,assignments)
		print(assignments)
		print(np.sum(assignments),len(assignments[assignments!=0]))
		i +=1
		t_max_mag = max(mags[assignments==0])
	return assignments

	
def one_mbc(lons,lats,deps,mags,times,ids,assignments):
	"""
	One iteration of magnitude based clustering

	assignements is an array that keeps track of the events 
	that have already been assigned
	"""
	# work only on unassigned earthquakes
	# find largest event
	ass = (assignments==0).astype(int)
	max_mag = np.max(mags*ass)
	max_arg = np.argmax(mags*ass)
	print('largest magnitude = ', max_mag)
	
	# estimate windows
	space_window, time_window = window(max_mag)

	t_lon = lons[max_arg]
	t_lat = lats[max_arg]
	t_dep = deps[max_arg]
	t_id  = ids[max_arg]
	t_time = times[max_arg]
	
	print(t_id, 'temporary Event ID')
	"""
	t_lons = lons[assignments==0]
	t_lats = lats[assignments==0]
	t_deps = deps[assignments==0]
	t_times = times[assignments==0]
	"""
	
	# calculate distances
	distances = calculate_distances(lats,lons,deps,t_lat,t_lon,t_dep)
	# calculate time delays
	delays    = calculate_delays(times,t_time)


	# assign new ids
	assignments[(assignments==0)&(distances<=space_window)&(delays<=time_window)&(delays>0)] = t_id
	# assign the mainshock 
	assignments[max_arg] = t_id	

	return assignments


def window(magnitude):
	"""
	function window: calculates the space-time window
	of aftershocks fro a given magnitude. This follows
	Gardner & Knoppoff, 1974 
	
	The parameters used come from running the code gardner&knopoff.py
	For the time window, these are the parameters from fitting
	a quadratic.
	Space window parameers come from a Weibull fit
	"""
	space_params = [9.73621447e+02, 1.62369152e+01, 8.46686624e-08, 8.93218484e+00]
	sp           = space_params
	time_params  = [ 1.75124875, -5.29020979, 22.62087912]
	tp           = time_params
	#name are flipped, solved below	

	time_window = sp[0] + (sp[1]-sp[0])*np.exp(-sp[2]*magnitude**sp[3])
	space_window  = tp[0]*magnitude**2 + tp[1]*magnitude + tp[2]

	print('space window = ',space_window, 'time window = ',time_window)
	return space_window, time_window 

def calculate_distances_old(lats,lons,deps,ref_lat,ref_lon,ref_dep):
	"""
	calculate the distaces between hypocenters, using the
	flat earth approximation
	"""
	R    = 6371 # km
	dlat = np.radians(lats - ref_lat)
	dlon = np.radians(lons - ref_lon)
	ddep = np.radians(deps - ref_dep)
	# calculate mean latitudes
	mlat = np.radians((lats + ref_lat) / 2)
	
	a         = np.sqrt(dlat**2 + (np.cos(mlat)*dlon)**2  )
	distances = np.sqrt((R*a)**2 + ddep**2)

	return distances

def calculate_distances(lats,lons,deps,ref_lat,ref_lon,ref_dep):
	"""
	Calculates the distances using first:
	a spherical to cartesian coordinate conversion
	and then the L2-norm
	"""
	R    = 6371 # km
	ref_dep = R - ref_dep 

	x_ref = ref_dep*np.sin(np.radians(90-ref_lat))*np.cos(np.radians(ref_lon))
	y_ref = ref_dep*np.sin(np.radians(90-ref_lat))*np.sin(np.radians(ref_lon))
	z_ref = ref_dep*np.cos(np.radians(90-ref_lat))
	
	x_ = (R-deps)*np.sin(np.radians(90-lats))*np.cos(np.radians(lons))
	y_ = (R-deps)*np.sin(np.radians(90-lats))*np.sin(np.radians(lons))
	z_ = (R-deps)*np.cos(np.radians(90-lats))

	distances = np.sqrt((x_ - x_ref)**2 + (y_ - y_ref)**2 + (z_ - z_ref)**2)

	return distances

def calculate_delays(times,ref_time):
	"""
	function calculate_delays
	take the inputs as strings, convert to UTCDateTime and calculate delays
	"""
	ut_times = [UT(time) for time in times]
	ref_ut_time = UT(ref_time)	

	delays = [ut_time - ref_ut_time for ut_time in ut_times]
	# delays are in seconds, convert to days
	delays = np.asarray(delays)
	delays = delays / 86400

	return delays
