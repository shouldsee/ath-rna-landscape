'''
Purpose:
   Visualise global transcriptome landscape for Arabidopsis thaliana.

Data from: 
- https://www.nature.com/articles/ncomms15309

'''

try:
	import scipy.io
except Exception as e:
	print(e)
from sklearn import manifold, datasets
from path import Path
import pickle
from spiper.types import Flow,Node, File, LinkFile
# from pand
import pandas as pd


from spiper.types import LoggedShellCommand

def download_tair10_defline(self,prefix,_output=['tsv']):
	CMD = ['curl','-sL','ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_functional_descriptions','>',self.output.tsv,]
	res = LoggedShellCommand(CMD)


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

def dump_columns(self,prefix, input=File, _output=['csv']):
	pd.DataFrame(_read_pandas(input).columns).to_csv(self.output.csv,index=0,header=None)


def sub_sample(self,prefix,
	input = File,
	_rowsep = 1,
	_colsep = 5,
	# rowsep=int,
	# colsep=int,
	_output=['pd_pk']):
	df = pd.read_pickle(input)
	df = df.iloc[::_rowsep, ::_colsep]
	df.to_pickle(self.output.pd_pk)

def sub_sample_genes(self,prefix,
	input = File,
	_rowsep = 5,
	_colsep = 1,
	# rowsep=int,
	# colsep=int,
	_output=['pd_pk']):
	df = _read_pandas(input)
	df = df.iloc[::_rowsep, ::_colsep]
	df.to_pickle(self.output.pd_pk)

def pandas_meannorm(self,prefix,input=File,_output=['pd_pk']):
	import numpy as np
	df = _read_pandas(input)
	df.loc[:] = df.values - np.mean(df.values,axis=1,keepdims=1)
	df.to_pickle(self.output.pd_pk)


def filter3(self,prefix, input_pd=File,data_pd=File, _output=['csv','pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[0] <0 )& (df.T[1] < 50 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_csv(self.output.csv)
	odf.to_pickle(self.output.pd_pk)

def filter5(self,prefix, input_pd=File,data_pd=File, _output=['csv','pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[1] - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_csv(self.output.csv)
	odf.to_pickle(self.output.pd_pk)

def filter6(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[0] < 300)
	# - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)

def filter7(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	# s = 'AT1G22770+AT1G01060'
	s = '''AT1G01060,AT3G09600,AT2G46830,AT4G15430,AT4G38960,AT3G02380,AT5G15850,AT5G17300,AT2G41250,AT1G64500,AT2G24540,AT3G12320,AT5G03555,AT5G52570,AT1G01520,AT2G31380,AT2G47490,AT5G18670,AT5G54130,AT1G69830,AT5G14760,AT1G10760,AT4G27360,AT1G32900,AT4G26670,AT5G15950,AT1G56300,AT4G16146,AT3G27170'''	
	# odf = _read_pandas(data_pd).loc[s.split('+')]
	odf = _read_pandas(data_pd).loc[s.strip(',').split(',')[:20]]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)

def filter_lhy(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	# s = 'AT1G22770+AT1G01060'
	s = '''AT1G01060,AT3G09600,AT2G46830,AT4G15430,AT4G38960,AT3G02380,AT5G15850,AT5G17300,AT2G41250,AT1G64500,AT2G24540,AT3G12320,AT5G03555,AT5G52570,AT1G01520,AT2G31380,AT2G47490,AT5G18670,AT5G54130,AT1G69830,AT5G14760,AT1G10760,AT4G27360,AT1G32900,AT4G26670,AT5G15950,AT1G56300,AT4G16146,AT3G27170'''	
	# odf = _read_pandas(data_pd).loc[s.split('+')]
	odf = _read_pandas(data_pd).loc[s.strip(',').split(',')[:20]]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)



def filter10(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	_ = '''
	curl http://localhost:8080/gene/filter4_genes/AT1G02340/csv | head -n101 | tail -n+2 | cut -d',' -f1 | tr '\n' ',' | tee temp
	'''
	# s = 'AT1G22770+AT1G01060'
	s = 'AT1G02340,AT1G74840,AT1G09570,AT4G13030,AT3G10740,AT3G48530,AT4G28040,AT4G28880,AT4G01120,AT1G27150,AT1G19870,AT1G64110,AT1G22640,AT1G53190,AT1G27670,AT3G62860,AT3G29035,AT3G61260,AT2G31980,AT1G21920,AT4G34970,AT4G27260,AT2G05540,AT1G68020,AT1G23050,AT3G60030,AT4G34000,AT4G31380,AT1G78700,AT5G45310,AT2G22680,AT3G11690,AT3G57410,AT3G10550,AT1G23870,AT3G10770,AT5G63190,AT4G23870,AT2G40970,AT3G61060,AT2G47485,AT3G62090,AT3G04350,AT1G34110,AT3G56930,AT1G12780,AT5G02020,AT1G75540,AT3G23640,AT3G04730,AT3G06210,AT2G45830,AT3G48260,AT4G37220,AT5G15330,AT5G23670,AT1G02860,AT1G48840,AT5G19740,AT4G26140,AT4G40060,AT3G54510,AT3G55770,AT4G14130,AT5G49710,AT5G02090,AT1G29760,AT2G45820,AT3G16800,AT5G51150,AT5G47720,AT5G63640,AT3G07880,AT4G34950,AT2G33990,AT4G33150,AT4G28650,AT2G28110,AT3G26890,AT2G26710,AT1G33055,AT3G15070,AT4G39675,AT3G51910,AT2G28200,AT4G21210,AT3G10250,AT2G46800,AT2G31960,AT3G56250,AT2G21820,AT1G75040,AT3G21260,AT3G05900,AT1G25230,AT3G01472,AT3G01470,AT1G55810,AT5G18650,AT4G16780'
	# s = 'AT4G16780,AT1G78700,AT1G21920,AT1G02340,AT4G01000,AT5G02200,AT1G09250,AT3G57410,AT3G10770,AT5G66030,AT5G48160,AT4G36730,AT2G41750,AT3G16857,AT1G03610,AT1G23050,AT1G04300,AT1G05890,AT1G74840,AT3G61260,AT1G12780,AT3G15290,AT3G06210,AT4G28260,AT4G32280,AT3G19860,AT4G35770,AT1G05840,AT1G36070,AT5G47180,'
	odf = _read_pandas(data_pd).loc[s.strip(',').split(',')[:15]]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)	


def filter_athb2(self,prefix, input_pd=str,data_pd=File, _output=['pd_pk']):
	_ = '''
	ACC=AT4G16780
	curl http://localhost:8080/gene/filter4_genes/${ACC}/csv | head -n101 | tail -n+2 | cut -d',' -f1 | tr '\n' ',' | tee temp
	'''
	s = 'AT4G16780,AT1G78700,AT1G21920,AT1G02340,AT4G01000,AT5G02200,AT1G09250,AT3G57410,AT3G10770,AT5G66030,AT5G48160,AT4G36730,AT2G41750,AT3G16857,AT1G03610,AT1G23050,AT1G04300,AT1G05890,AT1G74840,AT3G61260,AT1G12780,AT3G15290,AT3G06210,AT4G28260,AT4G32280,AT3G19860,AT4G35770,AT1G05840,AT1G36070,AT5G47180,AT4G19190,AT2G26900,AT2G26210,AT2G37220,AT4G28880,AT2G42890,AT3G04790,AT5G57660,AT1G06570,AT3G04560,AT4G10925,AT3G14750,AT2G26430,AT1G03140,AT5G45310,AT1G69220,AT5G13010,AT3G55770,AT3G22961,AT3G10030,AT1G06400,AT3G51880,AT3G20770,AT1G15140,AT2G01100,AT2G43060,AT1G79700,AT4G39910,AT5G39660,AT4G27450,AT3G23030,AT4G16670,AT1G53240,AT4G16110,AT1G29760,AT5G66110,AT5G66180,AT4G21450,AT4G00165,AT1G03475,AT5G18630,AT1G33050,AT3G01520,AT4G24230,AT5G46180,AT2G28200,AT1G09570,AT5G05750,AT4G09830,AT4G34740,AT5G19750,AT5G55280,AT5G44260,AT3G60300,AT4G27260,AT3G17850,AT1G10830,AT5G59070,AT1G28960,AT5G49910,AT3G47610,AT3G61060,AT3G03870,AT3G01472,AT3G01470,AT5G57940,AT4G15780,AT3G51110,AT3G16800,AT2G36900,'
	acc1 = s.strip(',').split(',')[:25]
	s = '''AT1G01060,AT3G09600,AT2G46830,AT4G15430,AT4G38960,AT3G02380,AT5G15850,AT5G17300,AT2G41250,AT1G64500,AT2G24540,AT3G12320,AT5G03555,AT5G52570,AT1G01520,AT2G31380,AT2G47490,AT5G18670,AT5G54130,AT1G69830,AT5G14760,AT1G10760,AT4G27360,AT1G32900,AT4G26670,AT5G15950,AT1G56300,AT4G16146,AT3G27170'''	
	acc2 = s.strip(',').split(',')[:25]
	odf = _read_pandas(data_pd).loc[acc1+acc2]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)



def filter11(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[1] < 7)
	# - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)


def filter12(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[1] > 7) & (df.T[1] < 17)
	# - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)

def filter13(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[1] > -5) & (df.T[1] < 5) & (df.T[0] < 10)
	# - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)

def filter14(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	# colsel = (df.T[1] > 10)
	colsel = (df.T[1] < 10) &(df.T[0]<20) & (df.T[1] > -10)
	# colsel = (df.T[0] > 20)
	# colsel = (df.T[1] < -7.5) | (df.T[1] > 11)
	odf = _read_pandas(data_pd).loc[:, colsel]
	odf.to_pickle(self.output.pd_pk)

def filter15(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	# colsel = (df.T[1] > 10)
	colsel = (df.T[1] < 0)
	 # &(df.T[0]<20)
	# colsel = (df.T[0] > 20)
	# colsel = (df.T[1] < -7.5) | (df.T[1] > 11)
	odf = _read_pandas(data_pd).loc[:, colsel]
	odf.to_pickle(self.output.pd_pk)


def filter16(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	# colsel = (df.T[1] > 10)
	# colsel = (df.T[1] > 7)
	colsel = (abs(df.T[2]) < 5) &( abs(df.T[1]) < 5)
	 # &(df.T[0]<20)
	# colsel = (df.T[0] > 20)
	# colsel = (df.T[1] < -7.5) | (df.T[1] > 11)
	odf = _read_pandas(data_pd).loc[:, colsel]
	odf.to_pickle(self.output.pd_pk)


def filter8(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[1] - (0.5*df.T[0] +10) < 0)
	# - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)


def filter9(self,prefix, input_pd=File,data_pd=File, _output=['pd_pk']):
	df = _read_pandas(input_pd)
	colsel =  (df.T[2] < 0.)
	 # - (0.5*df.T[0] +10) < 0)
	# colsel =  (df.T[1] - (0.5*df.T[0] +10) < 0)
	# - df.T[0] >-400 )
	odf = _read_pandas(data_pd).loc[:, colsel]
	# fn = self.prefix_named+'.filtered3.csv'
	odf.to_pickle(self.output.pd_pk)



def filter4(self,prefix, input_pd=File, data_pd=File, _output=['pd_pk']):	
	df = _read_pandas(input_pd)
	# df0 = df.copy()
	xs =  df.loc['AT1G22770']
	xs =  xs - xs.mean()
	print(xs<-2.)
	odf = _read_pandas(data_pd).loc[:,~(xs < -2.)]
	# odf = _read_pandas(data_pd)
	# .loc[:,xs < -2.]
	odf.to_pickle(self.output.pd_pk)

def filter_leaves(self,prefix, input_pd=File, data_pd=File, _output=['pd_pk']):	
	df = _read_pandas(input_pd)
	colsel = df['SAMPLE_ID'][(
		(
	 df['WORDS'].str.contains("leaf")     
	 |df['WORDS'].str.contains("leaves")   
	 |df['WORDS'].str.contains("leafs")    
	 |df['WORDS'].str.contains("seedling") 
	 |df['WORDS'].str.contains("seedlings")
	 )
	 &~df['WORDS'].str.contains("shoot")
	 &~df['WORDS'].str.contains("root")

	 &~df['WORDS'].str.contains("true dark whole seedling") ### this is a PIP-Seq, see https://www.ncbi.nlm.nih.gov/sra/?term=SRR1503852
	 &~df['WORDS'].str.contains("e-mtab-1499")### these are roots

	 &~df['WORDS'].str.contains("e-mtab-3279")### abnormal HOS1 AT2G39810 expression 
	 
	 &~df['WORDS'].str.contains("21 degrees australia 6dpi") ### abnormal HOS1 and LHY expression see https://www.ncbi.nlm.nih.gov/sra/?term=SRR1760178

	 &~df['WORDS'].str.contains("flowers") ###  these are flowers SRR360205
	 )]
	odf = _read_pandas(data_pd).loc[:,colsel]
	odf.to_pickle(self.output.pd_pk)

def filter_leaves_contaminated(self,prefix, input_pd=File, data_pd=File, _output=['pd_pk']):	
	df = _read_pandas(input_pd)
	colsel = df['SAMPLE_ID'][(
		(
	 df['WORDS'].str.contains("leaf")     
	 |df['WORDS'].str.contains("leaves")   
	 |df['WORDS'].str.contains("leafs")    
	 |df['WORDS'].str.contains("seedling") 
	 |df['WORDS'].str.contains("seedlings")
	 )
	 &~df['WORDS'].str.contains("shoot")
	 &~df['WORDS'].str.contains("root")

	 # &~df['WORDS'].str.contains("true dark whole seedling") ### this is a PIP-Seq, see https://www.ncbi.nlm.nih.gov/sra/?term=SRR1503852
	 # &~df['WORDS'].str.contains("e-mtab-1499")### these are roots

	 )]
	odf = _read_pandas(data_pd).loc[:,colsel]
	odf.to_pickle(self.output.pd_pk)



from spiper.types import DirtyKey
import collections
@Flow
def main(self, prefix, input =File,
	_output=[]):
	xml_prefix = self.prefix+'.xml/root'
	curr = self.runner(download_tair10_defline, prefix)
	curr = self.runner(LinkFile,   self.prefix_named +'.mat', input        )
	curr = self.runner(prep,         prefix,                    curr.prefix  )		
	_    = self.runner(dump_columns, prefix,                    curr.output.pd_pk )		
	curr = self.runner(sub_sample,   prefix,                    curr.output.pd_pk)
	curr = self.runner(logtpm_normalise,      prefix,          curr.output.pd_pk)
	_ = self.config_runner(tag='all')(logtpm_normalise,      prefix,          self.subflow['prep'].output.pd_pk)

	curr = self.runner(project_sample_isomap, prefix,          curr.output.pd_pk,                  1, 10)
	_    = self.config_runner(tag='project_sample_isomap')(fetch_sample_attr,       prefix,        xml_prefix,   
		self.subflow['project_sample_isomap'].output.pd_pk)
	_    = self.config_runner(tag='project_sample_isomap')(prepare_csv,             prefix, 			
		self.subflow['fetch_sample_attr-project_sample_isomap'].output.pd_pk,
		self.subflow['project_sample_isomap'].output.pd_pk,
		)	                  

	self.runner(filter_leaves,  prefix,  'root.fetch_sample_attr.pd_pk',   self.subflow['logtpm_normalise-all'].output.pd_pk)

	self.runner(filter3, prefix, curr.output.pd_pk,   self.subflow['logtpm_normalise'].output.pd_pk)
	self.runner(filter4, prefix, self.subflow['filter3'].output.pd_pk, self.subflow['filter3'].output.pd_pk )
	_    = self.runner(sub_sample_genes,      prefix,          self.subflow['filter3'].output.pd_pk)
	_    = self.config_runner(tag='all')\
					  (sub_sample_genes,      prefix,          self.subflow['logtpm_normalise'].output.pd_pk)
	_    = self.config_runner(tag='all')\
					  (pandas_meannorm,       prefix,          self.subflow['sub_sample_genes-all'].output.pd_pk)

	_    = self.config_runner(tag='all_genes')\
					  (pandas_meannorm,       prefix,          self.subflow['logtpm_normalise'].output.pd_pk)


	curr_genes_sample\
	     = self.config_runner(tag="filter3_genes")\
	     			  (pandas_meannorm,       prefix,          self.subflow['filter3'].output.pd_pk,)

	_    = self.config_runner(tag="filter4_genes")\
	     			  (pandas_meannorm,       prefix,          self.subflow['filter4'].output.pd_pk,)

	curr_genes_sample\
	     = self.config_runner(tag="filter3")\
	     			  (pandas_meannorm,       prefix,          self.subflow['sub_sample_genes'].output.pd_pk,)
	curr_genes = self.config_runner(tag='genes_meannorm_filter3')\
	                  (project_sample_isomap, prefix,          curr_genes_sample.output.pd_pk,     0, 10)


	curr_genes = self.config_runner(tag='filter4')\
	                  (project_sample_isomap, prefix,          self.subflow['filter4'].output.pd_pk,     1, 5)
	_ = self.runner(filter5, prefix, curr_genes.output.pd_pk,  self.subflow['filter4'].output.pd_pk)
	curr_genes = self.config_runner(tag='filter5')\
	                  (project_sample_isomap, prefix,          self.subflow['filter5'].output.pd_pk,     
	                  	1, 10)
	_ = self.runner(filter6, prefix, curr_genes.output.pd_pk,  self.subflow['filter5'].output.pd_pk)
	# _ = self.runner(filter7, prefix, curr_genes.output.pd_pk,  self.subflow['filter3'].output.pd_pk)
	_ = self.runner(filter7, prefix, curr_genes.output.pd_pk,  self.subflow['logtpm_normalise-all'].output.pd_pk)

	curr_genes = self.config_runner(tag='filter5_pca')\
	                  (project_sample_project, prefix,  'pca',        self.subflow['filter5'].output.pd_pk,     
	                  	1, 3,10)
	curr_genes = self.config_runner(tag='filter6')\
	                  (project_sample_isomap, prefix,          self.subflow['filter6'].output.pd_pk,     
	                  	1, 10)

	curr_genes = self.config_runner(tag='filter7')\
	                  (project_sample_isomap, prefix,          self.subflow['filter7'].output.pd_pk,     
	                  	1, 5)
	_ = self.runner(filter8, prefix, curr_genes.output.pd_pk,  self.subflow['filter7'].output.pd_pk)
	# _ = self.runner(filter8, prefix, curr_genes.output.pd_pk,  self.subflow['filter7'].output.pd_pk)

	curr_genes = self.config_runner(tag='filter7_pca')\
	                  (project_sample_project, prefix,  'pca',        self.subflow['filter7'].output.pd_pk,     
	                  	1, 3, 10)
	curr_genes = self.config_runner(tag='filter8_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter8'].output.pd_pk,     
	                  	1, 3, 5)
	_ = self.runner(filter9, prefix, curr_genes.output.pd_pk,  self.subflow['filter8'].output.pd_pk)

	curr_genes = self.config_runner(tag='filter8_pca')\
	                  (project_sample_project, prefix,   'pca',     self.subflow['filter8'].output.pd_pk,     
	                  	1, 3, 5)
	curr_genes = self.config_runner(tag='filter9_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter9'].output.pd_pk,     
	                  	1, 3, 5)
	curr_genes = self.config_runner(tag='filter9_pca')\
	                  (project_sample_project, prefix,   'pca',     self.subflow['filter9'].output.pd_pk,     
	                  	1, 3, 5)


	_ = self.runner(filter10, prefix, curr_genes.output.pd_pk,  self.subflow['filter3'].output.pd_pk)
	# _ = self.runner(filter10, prefix, curr_genes.output.pd_pk,  self.subflow['logtpm_normalise-all'].output.pd_pk)
	curr_genes = self.config_runner(tag='filter10_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter10'].output.pd_pk,     
	                  	1, 3, 5)

	_ = self.runner(filter11, prefix, curr_genes.output.pd_pk,  self.subflow['filter10'].output.pd_pk)
	curr_genes = self.config_runner(tag='filter11_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter11'].output.pd_pk,     
	                  	1, 3, 5)

	_ = self.runner(filter_athb2, prefix,  None, self.subflow['filter3'].output.pd_pk)
	# _ = self.runner(filter10, prefix, curr_genes.output.pd_pk,  self.subflow['logtpm_normalise-all'].output.pd_pk)
	curr_genes = self.config_runner(tag='filter_athb2_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter_athb2'].output.pd_pk,     
	                  	1, 3, 5)



	_    = self.config_runner(tag='filter_athb2')(fetch_sample_attr,       prefix,        xml_prefix,   
		self.subflow['filter_athb2'].output.pd_pk)
	curr = self.config_runner(tag='filter_athb2')(prepare_csv,             prefix, 			
		self.subflow['fetch_sample_attr-filter_athb2'].output.pd_pk,
		self.subflow['project_sample_project-filter_athb2_isomap'].output.pd_pk,
		)	                  



	# filter_lhy = filter7
	_ = self.runner(filter_lhy, prefix,  'NA',  self.subflow['filter3'].output.pd_pk)
	curr_genes = self.config_runner(tag='filter_lhy_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter_lhy'].output.pd_pk,     
	                  	1, 3, 5)
	_    = self.config_runner(tag='filter_lhy')(fetch_sample_attr,       prefix,        xml_prefix,   
		self.subflow['filter_lhy'].output.pd_pk)
	curr = self.config_runner(tag='filter_lhy')(prepare_csv,             prefix, 			
		self.subflow['fetch_sample_attr-filter_lhy'].output.pd_pk,
		self.subflow['project_sample_project-filter_lhy_isomap'].output.pd_pk,
		)


	# _ = self.runner(filter12, prefix, curr_genes.output.pd_pk,  self.subflow['filter10'].output.pd_pk)
	# curr_genes = self.config_runner(tag='filter12_isomap')\
	#                   (project_sample_project, prefix,   'isomap',     self.subflow['filter12'].output.pd_pk,     
	#                   	1, 3, 5)
	# _ = self.runner(filter13,  prefix,  curr_genes.output.pd_pk,  self.subflow['filter11'].output.pd_pk)
	# curr_genes = self.config_runner(tag='filter13_isomap')\
	#                   (project_sample_project, prefix,   'isomap',     self.subflow['filter13'].output.pd_pk,     
	#                   	1, 3, 5)

	_ = self.runner(filter14,  prefix,  
		self.subflow['project_sample_project-filter10_isomap'].output.pd_pk,  
		self.subflow['filter10'].output.pd_pk)		
	curr_genes = self.config_runner(tag='filter14_isomap')\
	                  (project_sample_project, prefix,   'isomap',     self.subflow['filter14'].output.pd_pk,     
	                  	1, 3, 5)

	# _    = self.config_runner(tag='filter14')(fetch_sample_attr,       prefix,        xml_prefix,   
	# 	self.subflow['filter14'].output.pd_pk)
	# curr = self.config_runner(tag='filter14')(prepare_csv,             prefix, 			
	# 	self.subflow['fetch_sample_attr-filter14'].output.pd_pk,
	# 	self.subflow['filter14'].output.pd_pk,
	# 	)

	# _ = self.runner(filter15,  prefix,  
	# 	self.subflow['project_sample_project-filter14_isomap'].output.pd_pk,  
	# 	self.subflow['filter14'].output.pd_pk)		
	# curr_genes = self.config_runner(tag='filter15_isomap')\
	#                   (project_sample_project, prefix,   'isomap',     self.subflow['filter15'].output.pd_pk,     
	#                   	1, 3, 5)
	# # self.config_runner(tag='filter16')(filter, prefix,)
	# _ = self.runner(filter16,  prefix,  
	# 	self.subflow['project_sample_project-filter15_isomap'].output.pd_pk,  
	# 	self.subflow['filter15'].output.pd_pk)		
	# _ = self.config_runner(tag='filter16_isomap')\
	#                   (project_sample_project, prefix,   'isomap',     self.subflow['filter16'].output.pd_pk,     
	#                   	1, 3, 5)

	# _    = self.config_runner(tag='filter16')(fetch_sample_attr,       prefix,        xml_prefix,   
	# 	self.subflow['filter16'].output.pd_pk)
	# curr = self.config_runner(tag='filter16')(prepare_csv,             prefix, 			
	# 	self.subflow['fetch_sample_attr-filter16'].output.pd_pk,
	# 	self.subflow['filter16'].output.pd_pk,
	# 	)

	# _ = self.runner(filter11, prefix, curr_genes.output.pd_pk,  self.subflow['logtpm_normalise-all'].output.pd_pk)
	# curr_genes = self.config_runner(tag='filter11_isomap')\
	#                   (project_sample_project, prefix,   'isomap',     self.subflow['filter11'].output.pd_pk,     
	#                   	1, 3, 5)

	curr_genes = self.config_runner(tag='genes_meannorm_all')\
	                  (project_sample_isomap, prefix,          self.subflow['pandas_meannorm-all'].output.pd_pk,     0, 10)

	curr_genes = self.config_runner(tag='genes')\
	                  (project_sample_isomap, prefix,           self.subflow['sub_sample_genes'].output.pd_pk,     0, 10)




	if not self.runner.is_meta_run:

		# df = _read_pandas(self.subflow['project_sample_isomap-genes'].output.pd_pk)
		# df = _read_pandas(self.subflow['project_sample_isomap-genes_meannorm_all'].output.pd_pk)
		df = _read_pandas(self.subflow['project_sample_isomap-genes_meannorm_filter3'].output.pd_pk)
		ref= _read_pandas(self.subflow['download_tair10_defline'].output.tsv)
		ref.index=ref.index.str.split('.',1).str.get(0)
		ref = ref.loc[~ref.index.duplicated()]

		import matplotlib.pyplot as plt
		fig,axs = plt.subplots(1,2)
		ax = axs[0]
		ax.scatter(df[0],df[1], s=0.5)
		ax.axis('tight')
		ax.set_xlim(-500,500)
		ax.set_ylim(-500,500)

		ax = axs[1]
		ax.scatter(df[2],df[1], s=0.5)
		ax.axis('tight')
		ax.set_xlim(-500,500)
		ax.set_ylim(-500,500)

		plt.gcf().savefig('temp.png')

		import numpy as np
		acc = 'AT5G61670'
		acc = 'AT5G62430'
		acc = 'AT1G66100'
		df = _read_pandas(self.subflow['pandas_meannorm-all'].output.pd_pk)
		ref= _read_pandas(self.subflow['download_tair10_defline'].output.tsv)
		ref.index=ref.index.str.split('.',1).str.get(0)
		ref = ref.loc[~ref.index.duplicated()]

		X = df.values * df.loc[[acc]].values
		covX =  pd.DataFrame({'cov':np.mean(X,axis=1,keepdims=0),'rmsd':np.sqrt(np.mean(np.square(df.values),axis=1))},df.index)
		covX['covsq']=covX['cov']**2
		covX['r']    =covX['cov'] /  covX['rmsd'] /np.sqrt( np.mean(np.square( df.loc[acc].values)))
		covX['rsq']  =covX['r'] ** 2
		idx  =  np.argsort(-covX['rsq'].values,axis=None)[:30]
		genes=  covX.iloc[idx]

		fig = plt.figure(figsize=[8,8])
		plt.scatter(covX['rmsd'],covX['cov'],s=2.,alpha=0.5)
		fig.savefig('temp2.png')
		with open('temp.html','w') as f:
			f.write('<img src="%s"></img>'%'temp2.png')
			df = genes
			pd.set_option("max_colwidth", -1)
			df = df.merge(ref,left_index=True,right_index=True)
			df.to_html(f)
			# f.flush()
			# print(ref.head(5))
			# ref.head(5).to_html(f)


	if 0:
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
		_    = self.config_runner(tag=DirtyKey(fn))(fetch_sample_attr, prefix, xml_prefix,             fn)

	if 0:
		# self.runner(filter3, prefix, curr.output.pd_pk, self.subflow['sub_sample'].output.pd_pk)
		fn = self.subflow['sub_sample'].output.pd_pk
		# curr = self.config_runner(tag='sub_sample')(project_sample_isomap,   prefix,        fn)
		_    = self.config_runner(tag='sub_sample')(fetch_sample_attr,       prefix,        xml_prefix,   fn)
		curr = self.config_runner(tag='sub_sample')(prepare_csv,             prefix, 			
			self.subflow['fetch_sample_attr-sub_sample'].output.csv,
			curr.output.pd_pk
			)

		if not self.runner.is_meta_run:
			import pickle
			with open(curr.output.count_pk, 'rb') as f:
				keywords = [ x[0] for x in pickle.load(f).most_common(60)]
				lst = []
		else:
			keywords = ['rosette']
		with open(self.prefix_named+'.sub_sample.html','w') as fo:
			for _keyword in keywords:
				curr_plot = self.config_runner(tag='sub_sample_%s'%DirtyKey(_keyword))(plot_keyword, prefix + '.plots/root', curr.output.csv, [_keyword])
				fo.write('\n<br><h2>keyword: %s</h2>'%_keyword)
				fo.write('\n<br><img src="%s"></img>'%curr_plot.output.png.relpath(self.prefix.dirname()))
				fo.flush()		

	if 0:
		self.runner(filter3, prefix, curr.output.pd_pk, self.subflow['sub_sample'].output.pd_pk)
		fn = self.subflow['filter3'].output.csv
		curr = self.config_runner(tag='filter3')(project_sample_isomap,   prefix,        fn,           1, 10)
		_    = self.config_runner(tag='filter3')(fetch_sample_attr,       prefix,        xml_prefix,   fn)
		curr = self.config_runner(tag='filter3')(prepare_csv, prefix, 			
			self.subflow['fetch_sample_attr-filter3'].output.csv,
			curr.output.pd_pk
			)

		if not self.runner.is_meta_run:
			import pickle
			with open(curr.output.count_pk, 'rb') as f:
				keywords = [ x[0] for x in pickle.load(f).most_common(60)]
				lst = []
		else:
			keywords = ['rosette']
		with open(self.prefix_named+'.filter3.html','w') as fo:
			for _keyword in keywords:
				curr_plot = self.config_runner(tag='filter3_%s'%DirtyKey(_keyword))(plot_keyword, prefix + '.plots/root', curr.output.csv, [_keyword])
				fo.write('\n<br><h2>keyword: %s</h2>'%_keyword)
				fo.write('\n<br><img src="%s"></img>'%curr_plot.output.png.relpath(self.prefix.dirname()))
				fo.flush()
		# curr = self.config_runner(tag='filter3_test')(plot_keyword, prefix, curr.output.csv, ['seedlings'])
	return self


def plot_keyword(self,prefix, data_pandas=File, keywords=list, _output=['png']):

	df = _read_pandas(data_pandas)
	import matplotlib.pyplot as plt

	fig, axs = plt.subplots(1,1,figsize=[8,8])
	ax = axs
	xs = df['0']
	ys = df['1']
	# cs = [False for x in df['WORDS']]
	cs = [any([_x in keywords for _x in str(x).split()]) for x in df['WORDS']]
	# cs = [x for x in df['WORDS'] any([_x in keywords for _x in x.split()])]
	plt.set_cmap('Set1')
	ax.scatter(xs,ys, c = cs, cmap=plt.cm.Set1,s= 5.)
	ax.axis('tight')
	ax.legend()	
	ax.set_title('keywords=%r'%keywords)
	fig.savefig(self.output.png)



def prepare_csv(self,prefix, meta_csv = File, data_csv = File, _output=['csv','count_pk']):
	if self.runner.is_meta_run:
		return 
	df1 = _read_pandas( data_csv).T
	df2 = _read_pandas( meta_csv).set_index('SAMPLE_ID')
	df = pd.concat([df1,df2[['source_name','tissue','WORDS']]],axis=1)
	df['WORDS'] = df['WORDS'].fillna('NA').astype(str)
	ct = collections.Counter()
	for x in df['WORDS']:
		# pprint(x)	
		ct.update(x.split())
	# pprint(ct.most_common(20))
	with open(self.output.count_pk,'wb') as f: pickle.dump(ct,f)
	df.to_csv( self.output.csv)

def logtpm_normalise(self,prefix, input=File, _output=['pd_pk']):
	import numpy as np
	df = _read_pandas(input)
	X = df.values
	X = np.log2(1+X)
	X = X - np.log2(np.sum(np.power(2,X),axis=0,keepdims=1))
	df.loc[:,:] = X
	df.to_pickle(self.output.pd_pk)



def project_sample_project(self, 
	prefix, 
	method =str, input_pk=File, transpose=int,n_components=int, n_neighbors=int,_output=['pd_pk','pd_csv','png']):
	from time import time
	import matplotlib.pyplot as plt
	# df = pd.read_pickle(input_pk)
	df = _read_pandas(input_pk)

	# n_neighbors = 10

	from collections import OrderedDict
	from functools import partial
	import numpy as np
	from sklearn import decomposition
	methods = OrderedDict()
	methods['isomap'] = manifold.Isomap(n_neighbors, n_components)
	methods['pca']    =decomposition.PCA(n_components)

	X = df.values
	if transpose:
		X = X.T
	n_samples = len(X)
	n_columns = len(X.T)

	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i,columns, %i neighbors"
	             % (n_samples,n_columns, n_neighbors), fontsize=14)


	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i columns, %i neighbors"
	             % (n_samples,n_columns, n_neighbors), fontsize=14)

	assert len(methods)!=0
	if 1:
		i = 0 
		label = method
		# method = methods[label]

		# for i, (label, method) in enumerate(methods.items()):
		t0 = time()
		Y = methods[method].fit_transform(X)
		print(Y.shape)
		t1 = time()
		print("%s: %.2g sec" % (label, t1 - t0))
		ax = fig.add_subplot(2, 5, 2 + i + (i > 3))
		ax.scatter(Y[:, 0], Y[:, 1], cmap=plt.cm.Spectral,s= 0.5)
		ax.set_title("%s (%.2g sec)" % (label, t1 - t0))

		if len(Y.T) >= 3:
			ax.axis('tight')
			ax = fig.add_subplot(2, 5, 2 + i + (i > 3) + 1)
			ax.scatter(Y[:, 2], Y[:, 1], cmap=plt.cm.Spectral,s= 0.5)
			ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
			ax.axis('tight')

	names = [ '%s%s'%(method,x) for x in range(len(Y.T))]
	if transpose:
		odf = pd.DataFrame(Y.T, names, df.columns)
	else:
		odf = pd.DataFrame(Y,   df.index, names)
		# range(len(Y.T)))
	odf.to_pickle(self.output.pd_pk)
	odf.to_csv(self.output.pd_csv)
	plt.gcf().savefig(self.output.png)
	# 'temp.png')
	print('[done]')

def project_sample_isomap(self,prefix, input_pk=File, transpose=int, n_neighbors=int,_output=['pd_pk','pd_csv','png']):
	from time import time
	import matplotlib.pyplot as plt
	# df = pd.read_pickle(input_pk)
	df = _read_pandas(input_pk)

	# n_neighbors = 10
	n_components = 3

	from collections import OrderedDict
	from functools import partial
	import numpy as np
	methods = OrderedDict()
	methods['isomap'] = manifold.Isomap(n_neighbors, n_components)

	X = df.values
	if transpose:
		X = X.T
	n_samples = len(X)

	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (n_samples, n_neighbors), fontsize=14)


	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (n_samples, n_neighbors), fontsize=14)

	assert len(methods)!=0
	for i, (label, method) in enumerate(methods.items()):
	    t0 = time()
	    Y = method.fit_transform(X)
	    print(Y.shape)
	    t1 = time()
	    print("%s: %.2g sec" % (label, t1 - t0))
	    ax = fig.add_subplot(2, 5, 2 + i + (i > 3))
	    ax.scatter(Y[:, 0], Y[:, 1], cmap=plt.cm.Spectral,s= 0.5)
	    ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
	    ax.axis('tight')
	    ax = fig.add_subplot(2, 5, 2 + i + (i > 3) + 1)
	    ax.scatter(Y[:, 2], Y[:, 1], cmap=plt.cm.Spectral,s= 0.5)
	    ax.set_title("%s (%.2g sec)" % (label, t1 - t0))
	    ax.axis('tight')

	if transpose:
		odf = pd.DataFrame(Y.T, range(len(Y.T)), df.columns)
	else:
		odf = pd.DataFrame(Y,   df.index, range(len(Y.T)))
	odf.to_pickle(self.output.pd_pk)
	odf.to_csv(self.output.pd_csv)
	plt.gcf().savefig(self.output.png)
	# 'temp.png')
	print('[done]')


def _read_pandas(input_pk,**kw):
	if input_pk.endswith('pk'):
		df = pd.read_pickle(input_pk,**kw)
	elif input_pk.endswith('csv'):
		df = pd.read_csv(input_pk,index_col=[0],**kw)
	elif input_pk.endswith('tsv'):
		df = pd.read_csv(input_pk,index_col=[0],sep='\t',**kw)
	else:
		assert 0, (input_pk,)
	return df


def fetch_sample_xml(self, prefix, sample_id=str,_output=['xml','json']):
		from collections import OrderedDict
		import io,requests,json
		import xml.etree.ElementTree as ET
		acc= sample_id
		print('[fetching]%s'%acc)
		resp = requests.get("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=exp&db=sra&term={acc}".format(**locals()))
		with open(self.output.xml,'wb') as f:
			f.write(resp.content)		
		root = ET.parse(io.BytesIO(resp.content)).getroot()
		# root = ET.parse(acc+'.xml').getroot()
		sample_attrs = root.findall('EXPERIMENT_PACKAGE/SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
		sample_attrs = OrderedDict([(x.findtext('TAG'),x.findtext('VALUE')) for x in sample_attrs])
		sample_attrs['SAMPLE_ATTRS_JSON'] = json.dumps(sample_attrs)
		sample_attrs['SAMPLE_ID'] = acc

		with open(self.output.json,'w') as f:
			json.dump(sample_attrs,f)
		# sample_attrs.

def fetch_sample_attr(self, prefix, xml_cache_prefix=str,  input_pk = File,_output=['csv','pd_pk']):
	from collections import OrderedDict
	from spiper.types import LoggedShellCommand
	import xml.etree.ElementTree as ET
	import requests,io
	import json
	# accs = "SRR1046851 SRR1046866".split()
	df = _read_pandas(input_pk)
	# df = pd.read_pickle(input_pk) 
	accs = df.columns
	df_proto = []
	for acc in accs:
		# runner = self.config_runner(tag=acc)
		try:
			curr = self.config_runner(tag=acc)(fetch_sample_xml, xml_cache_prefix, acc)
			if not self.runner.is_meta_run:
				sample_attrs = json.loads(open(curr.output.json,'r').read())
				sp = ' '.join(json.loads(sample_attrs['SAMPLE_ATTRS_JSON']).values()).lower().split()
				lst = []
				for x in sp:
					x = x.strip(',')
					if x in ['seedlings','seedling']:
						x = 'seedling'
					if x in ['rosettes','rosette']:
						x = 'rosette'
					if x in ['leaves', 'leaf', 'leafs']:
						x = 'leaf'
					lst.append(x)
				sample_attrs['WORDS'] = ' '.join(lst)
				# sample_attrs['WORDS'] = ' '.join([x for x in sample_attrs] )
				df_proto.append( pd.Series(sample_attrs))
		except Exception as e:
			print(e)
			import time
			# time.sleep(5)
	if not self.runner.is_meta_run:
		odf = pd.concat(df_proto,axis=1).T 
		# if df_proto else pd.DataFrame()
	else:
		odf = pd.DataFrame([],columns=['source_name','tissue','WORDS'])
	odf.__setitem__('source_name',odf.get('source_name',['NA']*len(odf)))
	odf.__setitem__('tissue',odf.get('tissue',['NA']*len(odf)))
	odf.__setitem__('WORDS',odf.get('WORDS',['']*len(odf)))
	odf.to_csv(self.output.csv)
	odf.to_pickle(self.output.pd_pk)




if __name__=='__main__':
	from spiper.runner import get_changed_files, cache_run
	from pprint import pprint
	tups = (main, '$PWD/root', 'NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered_public.mat')
	pprint(get_changed_files(*tups))

	runner = cache_run
	# runner = get_changed_files
	pprint(runner(*tups))
	print('[script_finished]')





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
	n_samples = len(X.T) 


	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (n_samples, n_neighbors), fontsize=14)


	fig = plt.figure(figsize=(15, 8))
	fig.suptitle("Manifold Learning with %i points, %i neighbors"
	             % (n_samples, n_neighbors), fontsize=14)



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
