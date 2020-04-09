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
		df = pd.read_pickle(input_pk)
	elif input_pk.endswith('csv'):
		df = pd.read_csv(input_pk,index_col=[0])
	elif input_pk.endswith('list'):
		df = pd.read_csv(input_pk,index_col=[0],header=None).index
	elif input_pk.endswith('tsv'):
		df = pd.read_csv(input_pk,index_col=[0],sep='\t')
	else:
		assert 0, (input_pk,)
	return df

def get_data_df(data_id):
	fn = 'root.{data_id}.pd_pk'.format(**locals())
	if Path(fn).isfile():
		df = _read_pandas(fn)
		return df
	fn = 'root.pandas_meannorm-{data_id}.pd_pk'.format(**locals())
	if Path(fn).isfile():
		df = _read_pandas(fn)
		return df
	else:
		assert 0,(fn,	)



# def get_dataset
def get_isomap(data_id, gene_id, meta_id=None):
	import base64
	from spiper import jinja2_format
	from sklearn.manifold import Isomap

	import time
	t0 = time.time()
	gene_id = GeneId(gene_id)
	n_neighbors = 10
	# xml_prefix = File('$PWD/root.xml/root').expand()
	root_prefix  = File('$PWD/root').expand()
	# DIR = 
	gene_id = GeneId(gene_id)

	# df = _show_gene_csv(gene_id, data_id)
	ref      = _read_pandas('root.download_tair10_defline.tsv')
	ref.index= ref.index.str.split('.',1).str.get(0)
	ref = ref.loc[~ref.index.duplicated()]
	_get_ref_df = ref

	df0 = df = get_data_df(data_id)

	with Path('tempdir.%s'%str(time.time())).makedirs_p().realpath() as cdir:
		####### filter by correlation
		### calculation
		X = df.values * df.loc[[gene_id]].values
		covX =  pd.DataFrame({'cov':np.mean(X,axis=1,keepdims=0),'rmsd':np.sqrt(np.mean(np.square(df.values),axis=1))},df.index)
		covX['r']    =covX['cov'].values /  covX['rmsd'].values /np.sqrt( np.mean(np.square( df.loc[gene_id].values)))

		idx  =  np.argsort(-abs(covX['r']).values,axis=None)
		genes=  covX.iloc[idx]
		df = genes
		df = pd.concat([df[['r','rmsd','cov']], _get_ref_df.reindex(df.index) ],axis=1)
		df.index.name = 'gene_id'

		### filtering
		accs= df.index[:25]
		n_feats = len(accs)
		df = df0
		df = df.loc[accs]

		prefix = cdir/'root'
		fn_pk = Path('temp.pk').realpath()
		df.to_pickle( fn_pk )


		### calculating isomap

		X = df.values.T
		n_components = 3
		Y = Isomap(n_neighbors=n_neighbors, n_components=n_components).fit_transform(X)
		df2 = pd.DataFrame(Y.T,['isomap%d'%i for i in range(n_components)], df.columns)

		# print(df2.head())

		# df2.to_pickle('df2.pk'); input_pk='df2.pk'

		# import metacsv_ath_rnaseq.utils
		# db = metacsv_ath_rnaseq.utils.dict_from_db_branch( "shouldsee", "metacsv-ath-rnaseq", "data")

		df_data = df2.T
		df_meta = _read_pandas(root_prefix+'.fetch_sample_attr.pd_pk').set_index("SAMPLE_ID")
		df_meta = df_meta[['source_name','tissue','WORDS']]

		df = pd.concat([df_data,df_meta],axis=1)
		df['WORDS'] = df['WORDS'].fillna('NA').astype(str)
		odf = df.reindex(df_data.index)
		myjson = odf.reset_index().to_json(orient='records')
		myTable = odf.to_html(table_id='myTable')

	cdir.rmtree();

	f = io.StringIO()
	pngBuffer = io.BytesIO()

	import matplotlib.pyplot as plt
	fig,axs = plt.subplots(1,2,figsize=[12,6])
	axs = axs.ravel()
	axs[0].scatter(odf['isomap0'], odf['isomap1'], s=4., )
	axs[1].scatter(odf['isomap2'], odf['isomap1'], s=4., )
	fig.suptitle('Isomap for %s points using %s features and %s neighbors'%(len(odf),n_feats,n_neighbors))
	fig.savefig(pngBuffer,format='png')
	pngBuffer.seek(0);
	pngString= base64.b64encode(pngBuffer.read()).decode('utf8')
	f.write('<h2>'+gene_id+'</h2>')
	f.write('<h2>'+"%.3fs"%(time.time() -t0) +'</h2>')
	# f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))
	with open('get_isomap.html','r') as ft:
		f.write(jinja2_format(ft.read(),**locals()))
	f.seek(0);
	return StreamingResponse(f,media_type='text/html')  