from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, Response, StreamingResponse
app = FastAPI()

app.mount("/static", StaticFiles(directory="."), name="static")

@app.get("/")
def read_root():
    return {"app": "ath_network"}

_ = '''
Genevestigator
Genemania
### showing necessity of computing correlation within a tissue
http://localhost:8080/gene_scatter/all_genes/AT5G61270+AT5G14740
http://localhost:8080/gene_scatter/all_genes/AT5G61270+AT2G24190
http://localhost:8080/gene_scatter/all_genes/AT1G22770+AT3G47470
'''

import io
import numpy as np
# from ath_network import _read_pandas
import pandas as pd
def _read_pandas(input_pk):
	if input_pk.endswith('pk'):
		df = pd.read_pickle(input_pk)
	elif input_pk.endswith('csv'):
		df = pd.read_csv(input_pk,index_col=[0])
	elif input_pk.endswith('tsv'):
		df = pd.read_csv(input_pk,index_col=[0],sep='\t')
	else:
		assert 0, (input_pk,)
	return df

# 	AT3G49530

import matplotlib.pyplot as plt
import base64


ref =_read_pandas('root.download_tair10_defline.tsv')
ref.index=ref.index.str.split('.',1).str.get(0)
ref = ref.loc[~ref.index.duplicated()]
_get_ref_df = ref
def get_ref_df():
	ref =_read_pandas('root.download_tair10_defline.tsv')
	# df = _read_pandas(self.subflow['pandas_meannorm-all'].output.pd_pk)
	# ref= _read_pandas(self.subflow['download_tair10_defline'].output.tsv)
	ref.index=ref.index.str.split('.',1).str.get(0)
	ref = ref.loc[~ref.index.duplicated()]
	return ref

# AT5G61270
from path import Path
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


def calc_r(gene_id, data_id):
	df = get_data_df(data_id)

	X = df.values * df.loc[[gene_id]].values
	covX =  pd.DataFrame({'cov':np.mean(X,axis=1,keepdims=0),'rmsd':np.sqrt(np.mean(np.square(df.values),axis=1))},df.index)
	covX['covsq']=covX['cov']**2
	covX['r']    =covX['cov'] /  covX['rmsd'] /np.sqrt( np.mean(np.square( df.loc[gene_id].values)))
	covX['rsq']  =covX['r'] ** 2
	return covX

# '/home/user/.local/bin/uvicorn --port 8080 main:app'
import numpy as np
@app.get("/data_scatter/{data_ids}/{gene_id}")
def show_data_scatter(gene_id:str,data_ids:str):

	data_id1,data_id2 = data_ids.split('+')
	cov1 = calc_r(gene_id,data_id1)
	cov2 = calc_r(gene_id,data_id2)
	ref  = get_ref_df()
	N = len(cov1)

	fig,axs = plt.subplots(1,1,figsize=[8,8])
	ax = axs
	xs, ys = abs(cov1['r']),abs(cov2['r'])
	r = np.corrcoef(xs,ys)[0,1]
	plt.scatter(xs,ys, s=4,alpha=0.5)
	ax.set_xlabel(data_id1)
	ax.set_ylabel(data_id2)

	xlim = plt.gca().get_xlim()
	ylim = plt.gca().get_ylim()
	plt.plot(xlim, ylim,'b--')
	plt.plot(xlim, [x-0.5 for x in ylim],'b--')
	plt.plot(xlim, [x+0.5 for x in ylim],'b--')
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	plt.title('N=%s,R=%.3f'%(N,r))

	pngBuffer = io.BytesIO()
	fig.savefig(pngBuffer)
	pngBuffer.seek(0);


	f = io.StringIO()
	f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))

	df = pd.DataFrame({'r_diff':abs(cov1['r'])-abs(cov2['r'])}).merge(ref,left_index=True,right_index=True)
	df.sort_values('r_diff').head(10).to_html(f)
	df.sort_values('r_diff').tail(10).to_html(f)

	f.seek(0);
	return HTMLResponse(f.read())

@app.get("/gene_scatter/{data_id}/{gene_ids:str}")
def show_gene_scatter(gene_ids:str,data_id):	
	gene_ids = gene_ids.split('+')

	# print((gene_ids,))
	gene_id1, gene_id2 = gene_ids
	df = get_data_df(data_id)
	ref= get_ref_df()
	# df = _read_pandas('root.{data_id}.pd_pk'.format(**locals()))
	# ref =_read_pandas('root.download_tair10_defline.tsv')
	# # df = _read_pandas(self.subflow['pandas_meannorm-all'].output.pd_pk)
	# # ref= _read_pandas(self.subflow['download_tair10_defline'].output.tsv)
	# ref.index=ref.index.str.split('.',1).str.get(0)
	# ref = ref.loc[~ref.index.duplicated()]

	# X = df.values * df.loc[[acc]].values
	# genes=  covX.iloc[idx]

	xs, ys = df.loc[gene_ids].values
	N = len(xs)
	r = np.corrcoef(xs,ys)[0,1]

	fig = plt.figure(figsize=[8,8])
	plt.scatter(xs,ys,s=4.,alpha=0.5)
	# plt.scatter(covX['rmsd'],covX['cov'],s=2.,alpha=0.5)
	plt.plot(plt.gca().get_xlim(), plt.gca().get_ylim(),'b--')
	# plt.title('N=%s'%len(df.values.T))
	plt.title('N=%s,R=%.3f'%(N,r))
	pngBuffer = io.BytesIO()
	fig.savefig(pngBuffer)
	pngBuffer.seek(0);
	# 'data:image/png;base64,'pngBuffer.read()

	f = io.StringIO()
	# with open('temp.html','w') as f:
	f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))
	# %'temp2.png')
	# df = genes
	# pd.set_option("max_colwidth", -1)
	# df = pd.DataFrame(genes,index=genes.index.values)
	# df = df[['r']].merge(ref,left_index=True,right_index=True)
	# df.to_html(f)

	f.seek(0);
	return HTMLResponse(f.read())


def _show_gene_csv(gene_id, data_id):
	gene_id = GeneId(gene_id)
	covX = calc_r(gene_id,data_id)
	idx  =  np.argsort(-abs(covX['r']).values,axis=None)
	genes=  covX.iloc[idx]
	df = genes
	df = df[['r','rmsd','cov']].merge(_get_ref_df,left_index=True,right_index=True)
	df.index.name = 'gene_id'
	return df

@app.get("/gene/{data_id}/{gene_id}/csv")
def show_gene_csv(gene_id, data_id):
	gene_id = GeneId(gene_id)
	return Response(
	    _show_gene_csv(gene_id, data_id).to_csv(),
	    # processed_file,
	    media_type='text/csv',
	    # headers={'Content-Disposition': 'attachment;filename="bb-etiketten-samengevoegd_'+ current_date +'.csv"'}
	)	
    # return _show_gene_csv(gene_id, data_id).to_csv()


# f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))

@app.get("/gallery")
# /{data_id}/{gene_id}")
def gallery():
	f = io.StringIO()
	for fn in Path('').glob('*'):
		if fn.endswith('png'):
			f.write('<br><h3>%s</h3>'%fn)
			f.write('<img src="/static/%s"></img>'%(fn))
		# base64.b64encode(pngBuffer.read()).decode('utf8')))
	f.seek(0);
	return HTMLResponse(f.read(),)


# @app.get("/static/{fn}")
# def static(fn):
# 	 # StreamingResponse(io.BytesIO
# 	return StreamingResponse(open(fn,'rb'))
# 	# .read())

@app.get("/gene/{data_id}/{gene_id}")
def show_gene(gene_id,data_id):
	gene_id = GeneId(gene_id)

	import matplotlib.pyplot as plt
	import base64
	acc = gene_id
	# acc = 'AT5G61670'
	# acc = 'AT5G62430'
	# acc = 'AT1G66100'
	# df = _read_pandas('root.pandas_meannorm-all.pd_pk')
	# df = _read_pandas('root.pandas_meannorm-{data_id}.pd_pk'.format(**locals()))
	# df = _read_pandas('root.pandas_meannorm.pd_pk')

	covX = _show_gene_csv(gene_id,data_id)

	fig = plt.figure(figsize=[8,8])
	# plt.scatter(covX['rmsd'],covX['cov'],s=2.,alpha=0.5)
	plt.scatter(covX['rmsd'],covX['rmsd'] * covX['r'],s=2.,alpha=0.5)
	plt.plot([0,2],[0,2],'b--')
	plt.plot([0,2],[0,-2],'b--')
	plt.title('N=%s'%len(covX))
	pngBuffer = io.BytesIO()
	fig.savefig(pngBuffer)
	pngBuffer.seek(0);
	# 'data:image/png;base64,'pngBuffer.read()

	f = io.StringIO()
	# with open('temp.html','w') as f:
	f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))
	# %'temp2.png')
	df = covX
	pd.set_option("max_colwidth", -1)
	f.write('<h2>Positive correlation</h2><br>')
	df.head(30).to_html(f)
	f.write('<h2>negative correlation</h2><br>')
	df.sort_values('r').head(30).to_html(f)
	# df.tail(30).to_html(f)

	f.seek(0);
	return HTMLResponse(f.read())
		# f.flush()
		# print(ref.head(5))
		# ref.head(5).to_html(f)
# @app.get("/")

from ath_network import project_sample_project,fetch_sample_attr,prepare_csv
from spiper.runner import cache_run,force_run,get_changed_files
from spiper.types import File
from spiper import jinja2_format


def GeneId(s):
	return s.upper()
'''
import main
from importlib import reload;
reload(main)
from main import  show_isomap, project_sample_project,calc_r
from main import _show_gene_csv


import time
t0 = time.time()
show_isomap(*"filter6/AT2G25930".split("/"))
print("[done]%.3f"%(time.time()-t0))

%load_ext line_profiler
%lprun -f show_isomap show_isomap(*"filter_leaves/AT2G25930".split("/"))
%lprun -f show_isomap show_isomap(*"filter6/AT2G25930".split("/"))

%lprun -f show_isomap show_isomap(*"filter_leaves/AT2G25930".split("/"))

%lprun -f calc_r show_isomap(*"filter_leaves/AT2G25930".split("/"))

%lprun -f project_sample_project show_isomap(*"filter_leaves/AT2G25930".split("/"))

%lprun -f show_isomap show_isomap(*"filter_leaves/AT2G25930".split("/"))
%lprun -f project_sample_project show_isomap(*"filter_leaves/AT2G25930".split("/"))

%lprun -f _show_gene_csv show_isomap(*"filter_leaves/AT2G25930".split("/"))

'''

from get_isomap import get_isomap
@app.get("/isomap/{data_id}/{gene_id}")
def show_isomap(data_id,gene_id,):
	f = get_isomap(data_id,gene_id,None)
	return StreamingResponse(f,media_type='text/html')

# 104.238.149.248
# @app.get("/isomap/{data_id}/{gene_id}")
def show_isomap(data_id,gene_id,):
	import time
	t0 = time.time()
	gene_id = GeneId(gene_id)
	n_neighbors = 10
	# xml_prefix = File('$PWD/root.xml/root').expand()
	root_prefix  = File('$PWD/root').expand()
	# DIR = 
	gene_id = GeneId(gene_id)

	# df = _show_gene_csv(gene_id, data_id)

	ref =_read_pandas('root.download_tair10_defline.tsv')
	ref.index=ref.index.str.split('.',1).str.get(0)
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
		_runner = cache_run
		curr_isomap = _runner(project_sample_project, 
			prefix,         
			'isomap',
			fn_pk, 
			1,
			3,
			n_neighbors )
			# 10)
		### concatenating csv
		input_pk = curr_isomap.output.pd_pk
		fetched_pd_pk = root_prefix+'.fetch_sample_attr.pd_pk'
		curr   = _runner(prepare_csv,
		             prefix, 			
		             fetched_pd_pk,
		             input_pk,
			)	   
		odf = _read_pandas(curr.output.csv).reindex(_read_pandas(input_pk).columns)
		myjson = odf.reset_index().to_json(orient='records')
		myTable = odf.to_html(table_id='myTable')
	cdir.rmtree();
		
	f = io.StringIO()
	pngBuffer = io.BytesIO()
	import matplotlib.pyplot as plt
	fig,axs = plt.subplots(1,2,figsize=[12,6])
	axs = axs.ravel()
	axs[0].scatter(odf['isomap0'],odf['isomap1'], s=4., )
	axs[1].scatter(odf['isomap2'],odf['isomap1'], s=4., )
	fig.suptitle('Isomap for %s points using %s features and %s neighbors'%(len(odf),n_feats,n_neighbors))
	fig.savefig(pngBuffer,format='png')
	pngBuffer.seek(0);
	f.write('<h2>'+gene_id+'</h2>')
	f.write('<h2>'+"%.3fs"%(time.time() -t0) +'</h2>')
	f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))
	f.write(jinja2_format('''

<head>

<head>
    <title>{{gene_id}}</title>
    <!--
    <script src="https://requirejs.org/docs/release/2.3.5/minified/require.js" type="text/javascript"></script>
<script type="text/javascript">
require.config({
    paths: {
        DT: '//cdn.datatables.net/1.10.19/js/jquery.dataTables.min',
    }
});
require(["DT"], function(DT) {
    $(document).ready( () => {
        // Turn existing table addressed by its id `myTable` into datatable
        $("#myTable").DataTable();
    })
});
</script>
-->
    <script src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-3.4.1.min.js" type="text/javascript"></script>
    <script src="http://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js" type="text/javascript"></script>
    <link rel="stylesheet" type="text/css" href="//cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css">

	<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>

    <!-- The jQuery library is a prerequisite for all jqSuite products -->
    <!-- # <script type="text/ecmascript" src="../../../js/jquery.min.js"></script>  -->
    <script type="text/ecmascript" src="http://www.guriddo.net/demo/js/trirand/i18n/grid.locale-en.js"></script>
    <script type="text/ecmascript" src="http://www.guriddo.net/demo/js/trirand/src/jquery.jqGrid.js"></script>

	<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css"> 
    <link rel="stylesheet" type="text/css" media="screen" href="http://www.guriddo.net/demo/css/trirand/ui.jqgrid-bootstrap.css" />


</head>


<body>

{{myTable}}


<script type="text/javascript"> 
$(document).ready( () => {
        // Turn existing table addressed by its id `myTable` into datatable
        $("#myTable").DataTable();
    })

 </script>

<!--	

<div style="margin-left:20px">
    <table id="jqGrid"></table>
    <div id="jqGridPager"></div>
</div>


<script type="text/javascript"> 
	

$(document).ready(function () {
	var myjson={{myjson}};
		$("#jqGrid").jqGrid({
		datatype: "local",
		data: myjson,
		styleUI : 'Bootstrap',		
		autowidth:true,
//		shrinkToFit: true,
		colModel:[
		{label: 'SAMPLE_ID', name:  'SAMPLE_ID', width:150,fixed:true},
		{label: 'isomap_c0', name:  '0', width:120,fixed:true, sorttype: 'number' },
		{label: 'isomap_c1', name:  '1', width:120,fixed:true, sorttype: 'number' },
		{label: 'isomap_c2', name:  '2', width:120,fixed:true, sorttype: 'number' },
		{label: 'source_name', name:  'source_name', width:null},
		{label: 'tissue', name:  'tissue', width:null},
		{label: 'WORDS', name:  'WORDS', width:900,fixed:true}
		],
//		datatype: "json",
//		url: 'data.json',		
//		 colModel: [
//			{ label: 'Category Name', name: 'CategoryName', width: 75 },
//			{ label: 'Product Name', name: 'ProductName', width: 90 },
//			{ label: 'Country', name: 'Country', width: 100 },
//			{ label: 'Price', name: 'Price', width: 80, sorttype: 'integer' },
//			// sorttype is used only if the data is loaded locally or loadonce is set to true
//			{ label: 'Quantity', name: 'Quantity', width: 80, sorttype: 'number' }                   
//		],
        editurl: 'clientArray',
		onSelectRow: editRow,
		viewrecords: true, // show the current page, data rang and total records on the toolbar
		width: 780,
		height: 800,
		rowNum: 500,
		loadonce: true, // this is just for the demo
		pager: "#jqGridPager"
	});

    var lastSelection;

    function editRow(id) {
        if (id && id !== lastSelection) {
            var grid = $("#jqGrid");
            grid.jqGrid('restoreRow',lastSelection);
            grid.jqGrid('editRow',id, {keys: true} );
            lastSelection = id;
        }
    }
   
});</script>
-->


</body>		
		''',**locals()))
	f.seek(0);
	return StreamingResponse(f,media_type='text/html')           




@app.get("/r_scatter/{data_id}/{gene_ids}")
def show_r_scatter(gene_ids,data_id):
	gene_ids = gene_ids.split('+')
	gene_id1, gene_id2 = gene_ids
	ref = get_ref_df()

	# X = df.values * df.loc[[acc]].values
	# genes=  covX.iloc[idx]
	covX1 = calc_r(gene_id1, data_id)
	covX2 = calc_r(gene_id2, data_id)

	xs,ys = abs(covX1['r']),abs(covX2['r'])
	N = len(xs)
	r = np.corrcoef(xs,ys)[0,1]

	fig,axs = plt.subplots(1,1, figsize=[8,8])
	ax = axs
	plt.scatter(xs,ys,s=4.,alpha=0.5)
	# plt.scatter(covX['rmsd'],covX['cov'],s=2.,alpha=0.5)
	plt.plot(plt.gca().get_xlim(), plt.gca().get_ylim(),'b--')
	ax.set_xlabel(gene_id1)
	ax.set_ylabel(gene_id2)
	# plt.title('N=%s'%len(df.values.T))
	plt.title('N=%s,R=%.3f'%(N,r))
	pngBuffer = io.BytesIO()
	fig.savefig(pngBuffer)
	pngBuffer.seek(0);
	# 'data:image/png;base64,'pngBuffer.read()

	f = io.StringIO()
	# with open('temp.html','w') as f:
	f.write('<img src="data:image/png;base64,%s"></img>'%(base64.b64encode(pngBuffer.read()).decode('utf8')))
	# %'temp2.png')
	# df = genes
	pd.set_option("max_colwidth", -1)
	# df = pd.DataFrame(genes,index=genes.index.values)
	# df = df[['r']].merge(ref,left_index=True,right_index=True)
	# df.to_html(f)
	df = pd.DataFrame({'r_diff': (xs -ys).round(2),'r_sum':(xs+ys).round(2),'xs':xs.round(2),'ys':ys.round(2)}).merge(ref,left_index=True,right_index=True)
	df.sort_values('r_diff').head(10).to_html(f)
	df.sort_values('r_diff').tail(10).to_html(f)
	df.sort_values('r_sum').tail(50).to_html(f)

	f.seek(0);
	return HTMLResponse(f.read())
