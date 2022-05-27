import os
import sys
import pandas as pd

# Adding location of source code to system path
# os.path.dirname(__file__) gives the directory of
# current file. Put in updated path if running script from elsewhere
# os.path joins all the folders of a path together in a
# system independent way (i.e. will work equally well on Windows, linux etc)
sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")


#os.getcwd() gets the path of where you are running from
outdir = os.path.join(os.getcwd(), "./output_test")
prefix = "test_new"

cscm = CICEROSCM(
    {
        "gaspam_file": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nyend": 2100,
        "nystart": 1750,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
     },
)
#        "sunvolc":1,
#        "rf_sun_data":pd.read_csv("/div/qbo/utrics/RadiativeForcing/RFforSCM/solar_erf_ar6.txt", header=None, skiprows=1, index_col=0)
#"emstart":2000,
#"rf_sun_file": "/div/qbo/utrics/RadiativeForcing/RFforSCM/solar_erf_ar6.txt"
#"rf_volc_file":"/div/qbo/utrics/RadiativeForcing/RFforSCM/volcanic_erf_ar6.txt",

cscm._run({"output_folder": outdir, "output_prefix": prefix}, make_plot=False)

