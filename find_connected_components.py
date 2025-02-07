
# read the nearest neighbor distances, and the catalog, find the connected components and compare to 
# Zaliapin & BenZion 2013, number of singles and number of families


import sys
import pandas as pd
import networkx as nx

fname      ='hauksson_nnd_m4.csv'
fname      ='1990_2000_ncedc_nnd.csv'

fname = sys.argv[1]

df         = pd.read_csv(fname)

threshold = 1e-6

background = df[df['nn_eta']>=threshold]
clustered  = df[df['nn_eta']<threshold]

print(df)

# create the graph
G = nx.Graph()
# add the nodes
G.add_nodes_from(df['EventID'].to_list())
print(G)
# add the edges only from the clustered set
# this will be a tuple (EventID,nn_id)
for i,row in clustered.iterrows():
	temp = (row['EventID'],row['nn_id'])
	G.add_edge(*temp)
print(G)

subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
print(len(subgraphs))
# get only the subgraphs with two or more nodes
families=[]
for subgraph in subgraphs:
	if subgraph.number_of_nodes()>1:
		families.append(subgraph)

print(len(families))
for family in families:
	print(family.number_of_nodes(),family.nodes())

# write  a file with the connected components
df['family']=0
#df['familiy_mainshock']

# loop over families, in the appropriate indices, write the number of the family
for i,family in enumerate(families):
	# select the appropriate events
	tdf=df[df['EventID'].isin(family.nodes())]
	indices=tdf.index
	for index in indices:
		df.loc[index,'assignments'] = i
#df.loc[df['A'] > 2, 'B'] = new_val

outname = fname.split('.')[0]+'_assignements.csv'
df.to_csv(outname)
print('saved ',outname)






