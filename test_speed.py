import time
t0 = time.time()

# import main
# from main import  show_isomap

from importlib import reload;
import get_isomap; reload(get_isomap)
from get_isomap import get_isomap as show_isomap


# , project_sample_project,calc_r
# from importlib import reload;
# from main import  show_isomap, project_sample_project,calc_r
# from main import _show_gene_csv
# %load_ext line_profiler
# %lprun -f show_isomap 

print("[done]%.3f"%(time.time()-t0))
#show_isomap(*"filter6/AT2G25930".split("/"))
show_isomap(*"filter_leaves/AT2G25930".split("/"))
print("[done]%.3f"%(time.time()-t0))
