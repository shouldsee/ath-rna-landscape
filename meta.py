 ### reference: https://www.ncbi.nlm.nih.gov/books/NBK25499 

from sklearn import manifold, datasets
from path import Path
import pickle
from spiper.types import Flow,Node, File, LinkFile
# from pand
import pandas as pd
from ath_network import _read_pandas
from spiper.types import LoggedShellCommand
import lxml.etree 
from models import LocalSample

def xml_tostring(self, encoding='utf8',pretty_print=True,**kw):
	if not isinstance(self, (lxml.etree._Element,lxml.etree._ElementTree)):
		self = self._xml_root
	else:
		pass
	return lxml.etree.tostring(self, pretty_print=pretty_print,encoding=encoding,**kw).decode(encoding)
def _pprint(self):
	print(xml_tostring(self,pretty_print=True))


from eutils._internal.xmlfacades.base import Base
class ESummaryResult(Base):
	_root_tag = 'eSummaryResult'
class EFetchResult(Base):
	_root_tag = 'eFetchResult'


from eutils._internal.client import Client, ESearchResult
import io,os
import warnings
def fetch_ncbi_sra_samples(index, ):
	'''
	Take a list of SRA PRIMARY_ID and returns xml output
	'''

	ec = Client(api_key=os.environ.get("NCBI_API_KEY", None))
	esr = ec._qs.esearch(dict(db='sra',
		rettype='xml',
		retmax = len(index),	
		term=' OR '.join(['%s [ACCN]'%x for x in  index])))
	esr = ESearchResult(esr)

	if len(esr.ids)!=len(index):
		warnings.warn('[WARN] id/acc length mismathcing %d/%d'%(len(esr.ids),len(index)))
	ids = esr.ids

	esr = ec._qs.efetch(dict(
		id=','.join(map(str,ids)),
		db='sra',
		retmax=len(ids),
		))
	return io.BytesIO(esr)
def fetch_AccList_as_Xml(self,prefix,input=File,_output=['xml']):
	index = _read_pandas(input,header=None).index
	esr   = fetch_ncbi_sra_samples(index)
	with open(self.output.xml,'wb') as f:
		f.write(esr.read())


@Flow
def fetch_AccList_as_SimpleCsv(self, prefix, input = File, _output=[
	'csv']):
	import os,io
	import warnings
	from eutils._internal.client import Client, ESearchResult
	import xmltodict
	import json
	if not self.runner.is_meta_run:
		print('[fetching] %s sra records'%len(list(open(input,'r').read().rstrip().splitlines())))
	index = _read_pandas(input,header=None).index
	curr = self.runner(fetch_AccList_as_Xml, self.prefix_named, input)
	if not self.runner.is_meta_run:
		jdata = xmltodict.parse(open(curr.output.xml,'rb'),force_list=lambda path,key,value: True)
		samples = jdata['EXPERIMENT_PACKAGE_SET'][0]['EXPERIMENT_PACKAGE']
		samples = [ LocalSample.from_ncbi_efetch(expt_package) for expt_package in samples]

		run_ids = sum([x.dict()['RUN_ID_LIST'] for x in samples],[])
		missing = [x for x  in index if x not in run_ids]
		assert len(missing)==0,misisng

		df = pd.concat([pd.Series(x.to_simple_dict()) for x in samples],axis=1).T
		df = df.drop_duplicates(subset=['SAMPLE_ID'])
		# .set_index('ID')
		df.to_csv(self.output.csv, index=0)
		print('[fetching] got %s records'%len(df) )
	return self


import lxml.etree as etree
@Flow
def main(self, prefix, csv_file = File, _output=[]):
	test_csv_file = fn = csv_file+'.test.csv'
	if 1:
		with open(csv_file,'r') as f:
			with open(fn,'w') as fo:
				fo.write('\n'.join(f.read().splitlines()[:100]))
	self.config_runner(tag='test')(fetch_AccList_as_SimpleCsv,       prefix, test_csv_file) 
	self.config_runner(tag='production')(fetch_AccList_as_SimpleCsv, prefix, csv_file) 
	return self


if __name__ == '__main__':	
	from spiper.runner import get_changed_files,cache_run
	tups = (main, '$PWD/root', 'root.dump_columns.csv')
	# 'root.dump_columns.csv')
	runner = get_changed_files
	runner(*tups)
	runner = cache_run
	runner(*tups)
	# for runner in [get_changed_files, ]