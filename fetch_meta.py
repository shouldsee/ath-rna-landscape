from ath_network import fetch_sample_attr
from spiper.types import Flow, File

@Flow
def main(self,prefix, input_pk=File, _output=[]):
	xml_prefix = self.prefix+'.xml/root'
	# fn = 'root.filter3.pk'
	_    = self.runner(fetch_sample_attr,       prefix,        xml_prefix,   
		input_pk)
		# fn)
		# self.subflow['project_sample_isomap'].output.pd_pk)
	return self

if __name__ == '__main__':
	from spiper.runner import get_changed_files, cache_run
	from pprint import pprint
	# tups = (main, '$PWD/root', 'root.filter3.pk')
	tups = (main, '$PWD/root', 'root.filter7.pd_pk')

	 # 'NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered_public.mat')
	pprint(get_changed_files(*tups))

	runner = cache_run
	pprint(runner(*tups))
	print('[script_finished]')

