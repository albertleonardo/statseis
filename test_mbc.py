
import mbc
import numpy as np
import pandas as pd


df = pd.read_csv('ncedc.csv')
print(df)


lons=df['Longitude'].to_numpy()
lats=df['Latitude'].to_numpy()
deps=df['Depth'].to_numpy()
mags=df['Magnitude'].to_numpy()
times=df['DateTime'].to_list()
ids=df['EventID'].to_numpy()


assignments=mbc.mbc(lons,lats,deps,mags,times,ids,min_mag=3)
assignemts=assignments.astype(int)

df['assignments']=assignments
df.to_csv('assigned.csv')



