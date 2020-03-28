'''
Purpose:
   Visualise global transcriptome landscape for Arabidopsis thaliana.

Data from: 
- https://www.nature.com/articles/ncomms15309

'''

import scipy.io
from sklearn import manifold, datasets
from path import Path
import pickle
from spiper.types import Flow,Node, File, LinkFile
# from pand
import pandas as pd



def prep(self, prefix, 
	input = File,
	_output = ['pd_pk','mat']
	):

	# '''
	# ### gdown https://drive.google.com/uc?id=0B252lj6tAx8XdjJISV9QWXJ5dTg
	# ###
	# ''' 
	# [todo:spiper] update code comments as appropriate but don't count as identity

	# fn = "NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered_public.mat"
	# fn = Path(fn)
	# pkf = (fn+'.pk')
	# if not pkf.isfile():
	mat = scipy.io.loadmat(input)
	sids = [x[0] for x in mat['sids'][0][:]]
	tids = [x[0][0] for x in mat['tids'][:]]
	df = pd.DataFrame(mat['Y'], tids, sids)
	df.to_pickle(self.output.pd_pk)
	# with open( self.output.pk,'wb') as f: 
	# 	pickle.dump( df,f)

	# else:
	# 	with open(pkf,'rb') as f:
	# 		mat = pickle.load(f)
	# X = mat['Y']

def sub_sample(self,prefix,
	input = File,
	_rowsep = 2,
	_colsep = 5,
	# rowsep=int,
	# colsep=int,
	_output=['pd_pk']):
	df = pd.read_pickle(input)
	df = df.iloc[::_rowsep, ::_colsep]
	df.to_pickle(self.output.pd_pk)



from spiper.types import DirtyKey
@Flow
def main(self, prefix, input =File,
	_output=[]):
	curr = self.runner(LinkFile,   self.prefix_named +'.mat', input        )
	curr = self.runner(prep,       prefix,                    curr.prefix  )		
	curr = self.runner(sub_sample, prefix,                    curr.output.pd_pk)
	curr = self.runner(project_sample_isomap,prefix,          curr.output.pd_pk)


	if 1:
	# if self.mock == 0:
		df = pd.read_pickle(curr.output.pd_pk)
		# odf = df.loc[:, (df.T[0] > 250 )& (df.T[1] < 0 )]
		# print(odf.shape)
		# odf.to_csv(self.prefix_named+'.filtered.csv')
		# fn = self.prefix_named+'.filtered.csv'
		# _    = self.runner(fetch_sample_attr, prefix,             fn)


		odf = df.loc[:, (df.T[0] <0 )& (df.T[1] < 0 )].T.head(50).T
		fn = self.prefix_named+'.filtered2.csv'
		odf.to_csv(fn)
		_    = self.config_runner(tag=DirtyKey(fn))(fetch_sample_attr, prefix,             fn)

		# _    = self.runner(fetch_sample_attr, prefix,             curr.output.pd_pk)
		# curr = self.runner
		# curr = self.runner(sub_sample, prefix, )
	return self



def project_sample_isomap(self,prefix, input_pk=File, _output=['pd_pk','pd_csv','png']):
	from time import time
	import matplotlib.pyplot as plt
	df = pd.read_pickle(input_pk)

	n_neighbors = 10
	n_components = 3

	from collections import OrderedDict
	from functools import partial
	import numpy as np
	methods = OrderedDict()
	methods['isomap'] = manifold.Isomap(n_neighbors, n_components)
	X = df.values
	X = np.log2(1+X)
	X = X - np.log2(np.sum(np.power(2,X),axis=0,keepdims=1))
	# X = X.T

	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (len(X), n_neighbors), fontsize=14)


	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (len(X), n_neighbors), fontsize=14)

	assert len(methods)!=0
	for i, (label, method) in enumerate(methods.items()):
	    t0 = time()
	    Y = method.fit_transform(X.T)
	    print(Y.shape)
	    t1 = time()
	    print("%s: %.2g sec" % (label, t1 - t0))
	    ax = fig.add_subplot(2, 5, 2 + i + (i > 3))
	    ax.scatter(Y[:, 0], Y[:, 1], cmap=plt.cm.Spectral,s= 0.3	,)
	    ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
	    ax.axis('tight')

	odf = pd.DataFrame(Y.T, range(len(Y.T)), df.columns)
	odf.to_pickle(self.output.pd_pk)
	odf.to_csv(self.output.pd_csv)
	plt.gcf().savefig(self.output.png)
	# 'temp.png')
	print('[done]')



def fetch_sample_attr(self, prefix,   input_pk = File,_output=['csv']):
	from spiper.types import LoggedShellCommand
	import xml.etree.ElementTree as ET
	# accs = "SRR1046851 SRR1046866".split()
	if input_pk.endswith('pk'):
		df = pd.read_pickle(input_pk)
	elif input_pk.endswith('csv'):
		df = pd.read_csv(input_pk,index_col=[0])
	else:
		assert 0, (input_pk,)
	# df = 
	# df = pd.read_pickle(input_pk) 
	accs = df.columns
	lst = []
	for acc in accs:
		print('[fetching]%s'%acc)
	# for acc in df.columns:
		res = LoggedShellCommand([
		'curl "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=exp&db=sra&term={acc}" > {acc}.xml'.format(**locals())
		])
		root = ET.parse(acc+'.xml').getroot()
		sample_attrs = root.findall('EXPERIMENT_PACKAGE/SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
		sample_attrs = {x.findtext('TAG'):x.findtext('VALUE') for x in sample_attrs}
		sample_attrs['SAMPLE_ID'] = acc
		lst.append( pd.Series(sample_attrs))
	odf = pd.concat(lst,axis=1).T
	odf.to_csv(self.output.csv)



if __name__=='__main__':
	from spiper.runner import get_changed_files, cache_run
	from pprint import pprint
	tups = (main, '$PWD/root', 'NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered_public.mat')
	# pprint(get_changed_files(*tups))
	runner = cache_run
	# runner = get_changed_files
	pprint(runner(*tups))





def old():


	from time import time
	import matplotlib.pyplot as plt
	df = pd.read_pickle(input_pk)

	n_neighbors = 10
	n_components = 3

	from collections import OrderedDict
	from functools import partial
	methods = OrderedDict()
	LLE = partial(manifold.LocallyLinearEmbedding,
	              n_neighbors, n_components, eigen_solver='auto')
	# methods['LLE'] = LLE(method='standard')
	methods['LTSA'] = LLE(method='ltsa')
	# methods['Hessian LLE'] = LLE(method='hessian')
	# methods['Modified LLE'] = LLE(method='modified')
	methods['Isomap'] = manifold.Isomap(n_neighbors, n_components)
	# methods['MDS'] = manifold.MDS(n_components, max_iter=100, n_init=1)
	# methods['SE'] = manifold.SpectralEmbedding(n_components=n_components,
	#                                            n_neighbors=n_neighbors)
	# methods['t-SNE'] = manifold.TSNE(n_components=n_components, init='pca',perplexity=80,
	#                                  random_state=0)

	import numpy as np

	X = np.log2(1+X)
	X = X - np.log2(np.sum(np.power(2,X),axis=0,keepdims=1))
	X = X[:28000:2,:3000:5]
	# X = X - np.mean(X,axis=1,keepdims=1)
	# X = X - np.mean(X,axis=0,keepdims=1)
	X = X.T

	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (len(X), n_neighbors), fontsize=14)


	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (len(X), n_neighbors), fontsize=14)



	for i, (label, method) in enumerate(methods.items()):
	    t0 = time()
	    Y = method.fit_transform(X)
	    print(Y.shape)
	    t1 = time()
	    print("%s: %.2g sec" % (label, t1 - t0))
	    ax = fig.add_subplot(2, 5, 2 + i + (i > 3))
	    ax.scatter(Y[:, 0], Y[:, 1], cmap=plt.cm.Spectral,s= 0.3	,)
	    ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
	    ax.axis('tight')

	plt.gcf().savefig('temp.png')

	print('[done]')
	# import pdb;pdb.set_trace();
