import pandas as pd
import sys

pd.read_pickle(sys.argv[1]).to_csv(sys.stdout)