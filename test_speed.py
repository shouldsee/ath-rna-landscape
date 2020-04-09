import main
from importlib import reload;
reload(main)
from main import  show_isomap, project_sample_project,calc_r
from main import _show_gene_csv


import time
t0 = time.time()
show_isomap(*"filter6/AT2G25930".split("/"))
print("[done]%.3f"%(time.time()-t0))
