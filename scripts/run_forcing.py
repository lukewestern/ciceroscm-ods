import os
import shutil
import sys
import cProfile, pstats, io
import pandas as pd
#from pstats import SortKey

sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))
from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
outdir = os.path.join(os.getcwd(), "output_test_forcing")
#pr = cProfile.Profile()
#pr.enable()

prefix = "test_lambda"

cscm = CICEROSCM(
    {
        "sunvolc": 0,
        "nyend": 2100,
        "forc_file": os.path.join(data_dir, "CO2_1pros.txt"),
    },
)


pamset_udm_default={
    "rlamdo": 16.0,
    "akapa": 0.634,
    "cpi": 0.4,
    "W": 4.0,
    "beto": 3.5,
    "threstemp": 7.0,
    "lambda": 0.540,
    "mixed": 60.0,
    "foan": 0.61,
    "foas": 0.81,
    "ebbeta": 0.0,
    "fnso": 0.7531,}

pamset_udm={
    "rlamdo": 16.0,
    "akapa": 0.634,
    "cpi": 0.4,
    "W": 4.0,
    "beto": 3.5,
    "threstemp": 7.0,
    "lambda": 0.840,
    "mixed": 60.0,
    "foan": 0.61,
    "foas": 0.81,
    "ebbeta": 0.0,
    "fnso": 0.7531,}


cscm._run({"output_folder": outdir,
           "output_prefix": prefix},
          pamset_udm =pamset_udm
)



#pr.disable()
#s = io.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
#ps.print_stats()
#print(s.getvalue())
