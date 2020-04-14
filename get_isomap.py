from path import Path
File = Path
def GeneId(s):
	return s.upper()

def StreamingResponse(*x,**kw):
	return x[0]
import io
import numpy as np
# from ath_network import _read_pandas
import pandas as pd
def _read_pandas(input_pk):
	if input_pk.endswith('pk'):
		with open(input_pk,'rb') as f:
			df = pd.read_pickle(f)
		# df = pd.read_pickle(input_pk)
	elif input_pk.endswith('csv'):
		df = pd.read_csv(input_pk,index_col=[0])
	elif input_pk.endswith('list'):
		df = pd.read_csv(input_pk,index_col=[0],header=None)
	elif input_pk.endswith('tsv'):
		df = pd.read_csv(input_pk,index_col=[0],sep='\t')
	else:
		assert 0, (input_pk,)
	return df

# def get_data_df(data_id):
# 	fn = 'root.{data_id}.pd_pk'.format(**locals())
# 	if Path(fn).isfile():
# 		df = _read_pandas(fn)
# 		return df
# 	fn = 'root.pandas_meannorm-{data_id}.pd_pk'.format(**locals())
# 	if Path(fn).isfile():
# 		df = _read_pandas(fn)
# 		return df
# 	else:
# 		assert 0,(fn,	)






import matplotlib.pyplot as plt
# def get_dataset
import base64
# from spiper import jinja2_format
from sklearn.manifold import Isomap, TSNE
from jinja2 import StrictUndefined, Template
from sklearn import manifold

def jinja2_format(s,**context):
	# d = context.copy()
	# d = __builtins__.copy()
	# d.update(context)
	# .update(__builtins__)
	return Template(s,undefined=StrictUndefined).render(**context)

def get_isomap(data_df, gene_id, meta_id=None):
	print('[getting_isomap]')
	import time
	t0 = time.time()
	gene_id = GeneId(gene_id)
	# xml_prefix = File('$PWD/root.xml/root').expand()
	root_prefix  = File('$PWD/root').expand()
	# DIR = 
	gene_id = GeneId(gene_id)

	# df = _show_gene_csv(gene_id, data_id)
	ref      = _read_pandas('root.download_tair10_defline.tsv')
	ref.index= ref.index.str.split('.',1).str.get(0)
	ref = ref.loc[~ref.index.duplicated()]
	_get_ref_df = ref

	# df0 = df = get_data_df(data_id)
	df0 = df = data_df
	# df  = df0.copy()
	df.loc[:] = df.values - df.values.mean(axis=1,keepdims=True)
	# M = df.values.mean( axis=1, keepdims=True)

	# df0 = df
	# df0 = df = df.astype("float32")

	with Path('tempdir.%s'%str(time.time())).makedirs_p().realpath() as cdir:
		####### filter by correlation
		### calculation
		X = df.values * df.loc[[gene_id]].values
		covX =  pd.DataFrame({
			'cov': np.mean(X,axis=1,keepdims=0),'rmsd':np.sqrt(np.mean(np.square(df.values),axis=1))},df.index)
		covX = covX.astype("float32")
		covX['r'] = covX['cov'].values/  covX['rmsd'].values /np.sqrt( np.mean(np.square( df.loc[gene_id].values)))
		# plt.hist(covX['r'],bins=50)
		fig = plt.figure()
		plt.hist(covX['r'],bins=pd.np.linspace(-1,1,50))
		plt.ylim(0,50)
		# plt.scatter(df.values[:,0],df.values[:,1])

		fig.savefig('../temp.png')
		print([abs(covX['r']).min(),abs(covX['r']).max()])
		idx  =  np.argsort(-abs(covX['r']).values,axis=None)
		genes=  covX.iloc[idx]
		df = genes
		df = pd.concat([df[['r','rmsd','cov']], _get_ref_df.reindex(df.index) ],axis=1)
		df.index.name = 'gene_id'

		### filtering
		# accs= df.index[:35]
		# accs= df.index[:25]
		# accs= df.index[:55]
		n_neighbors = 10
		n_components = 2
		accs= df.index[:35]
		# accs= df.index[:15]
		n_feats = len(accs)

		df = df0
		# del df,df0		
		# with (cdir.dirname()):	
		# 	df = get_data_df(data_id)
		df = df.loc[accs]
		### calculating isomap

		X = df.values.T
		'''
		Use Isomap in favour for TSNE.
		reason see bottom of https://leovan.me/cn/2018/03/manifold-learning/
		https://d33wubrfki0l68.cloudfront.net/ab8f4f486b9be87e87a25672564607ad72e46904/2c261/images/cn/2018-03-16-manifold-learning/severed-sphere.png
		'''
		# random_state = None
		# plotRadius = 30;
		#### hessian-lle, not stable
		# model =  manifold.LocallyLinearEmbedding(n_neighbors, n_components,
  #                                                  eigen_solver='dense',
  #                                                  method='hessian',
  #                                                  random_state=random_state)

		plotRadius = 30;
		model = manifold.Isomap(n_neighbors=n_neighbors, n_components=n_components)
		# plotRadius = 70;
		# Y = TSNE(n_components=n_components,init='pca').fit_transform(X)

		Y = model.fit_transform(X)

		df2 = pd.DataFrame(Y.T,['isomap%d'%i for i in range(n_components)], df.columns)

		# print(df2.head())

		# df2.to_pickle('df2.pk'); input_pk='df2.pk'

		# import metacsv_ath_rnaseq.utils
		# db = metacsv_ath_rnaseq.utils.dict_from_db_branch( "shouldsee", "metacsv-ath-rnaseq", "data")

		df_data = df2.T
		df_meta = _read_pandas(root_prefix+'.fetch_sample_attr.pd_pk').set_index("SAMPLE_ID")
		df_meta["EXPT_ID"] = _read_pandas(root_prefix.dirname()/"bq_example.py.csv").set_index("acc")["sra_study"].reindex(df_meta.index)
		df_meta = df_meta[["EXPT_ID",'source_name','tissue','WORDS']]

		df = pd.concat([df_data,df_meta],axis=1)
		df['WORDS'] = df['WORDS'].fillna('NA').astype(str)
		odf = df.reindex(df_data.index)
		myjson = odf.reset_index().to_json(orient='records')
		myTable = odf.to_html(table_id='myTable')

	cdir.rmtree();
	
	f = io.StringIO()
	with  io.BytesIO() as pngBuffer:

		fig,axs = plt.subplots(1,2,figsize=[12,6])
		axs = axs.ravel()
		axs[0].scatter(odf['isomap0'], odf['isomap1'], s=4., )
		# axs[1].scatter(odf['isomap2'], odf['isomap1'], s=4., )
		fig.suptitle('Isomap for %s points using %s features and %s neighbors'%(len(odf),n_feats,n_neighbors))
		fig.savefig(pngBuffer,format='png')
		pngBuffer.seek(0);
		pngString= base64.b64encode(pngBuffer.read()).decode('utf8')
		f.write('<h2>'+gene_id+'</h2>')
		f.write('<h2>'+"%.3fs"%(time.time() -t0) +'</h2>')
		f.write("%r<br>"%([abs(covX['r']).values.min(), abs(covX['r']).values.max()]))
		f.write("%r"%([sorted(-abs(covX['r']).values)[::-1][:10],sorted(-abs(covX['r']).values)[:10]]))
		# f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))
		with open('get_isomap.html','r') as ft:
			f.write(jinja2_format(ft.read(),**locals()))
		f.seek(0);
	for k in locals().keys():
		if k!='f':
			del locals()[k] 

	return f
		# return StreamingResponse(f.read(),media_type='text/html')  