import sys
import re
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import pandas.testing as pdt
import warnings
import xarray as xr
from scipy.stats import qmc, norm
try:
    from pandas.core.common import SettingWithCopyWarning
except:
    from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
sys.path.insert(0,os.path.join(os.getcwd(), '../', 'src'))
from ciceroscm import CICEROSCM
from ciceroscm.input_handler import read_inputfile,read_components,read_natural_emissions

nsamp = 1000
fileindicator = sys.argv[1]
gcase = "halons"

# These are all the gases that are considered. Can add more if available or needed.
groups = {"all":[],
          "CFCs": ["CFC-11", "CFC-12", "CFC-113", "CFC-114", "CFC-115","CFC-13"],
          "HCFCs": ["HCFC-22", "HCFC-141b", "HCFC-142b", "HCFC-124", "HCFC-133a", "HCFC-123"],
          "HFCs": ["HFC-23", "HFC-32", "HFC-134a", "HFC-143a", "HFC-125",
                   "HFC-152a", "HFC-4310mee", "HFC-227ea", "HFC-365mfc",
                   "HFC-236fa", "HFC-245fa"], 
          "Solvents": ["CCl4", "CH3CCl3"],
          "PFCs":["CF4", "C2F6", "C3F8", "C4F8", "C4F10","C6F14"],
          "Other":["SF6", "NF3","SO2F2"],
          "Halons": ["H-1301", "H-2402", "H-1211"],
          "CH3Br": ["CH3Br"]}

# Gases to ger ERF for
csmsp = groups["Halons"]

def to_df(cscm):
    """Convert the results of a CICEROSCM object to a DataFrame"""
    out=pd.concat([pd.DataFrame(v) for k, v in cscm.results.items()], axis = 1, keys = list(cscm.results.keys()))
    return out

def calc_beta(mol_mass):
    """Calculate the beta value for a given molecular mass"""
    return 1.7758620689655172e8 * mol_mass * 1e-9


def check_hfc(species):
    """Check if the species is an HFC and remove the hyphen"""
    if "HFC-" in species:
        sp_in = species.replace("-", "")
    else:
        sp_in = species
    return sp_in

test_data_dir = "/user/home/lw13938/work/ciceroscm/tests/test-data/"

# Read gas parameters
# gaspam =read_components(test_data_dir + '/gases_v1RCMIP.txt')
gaspam = read_components("/user/home/lw13938/work/ciceroscm/inputfiles/gases_WMO2022.txt")
df_nat_ch4 =read_natural_emissions(test_data_dir + '/natemis_ch4.txt','CH4')
df_nat_n2o =read_natural_emissions(test_data_dir + '/natemis_n2o.txt','N2O')

# Add missing emissions and update with AGAGE numbers where available
df_ssp2_conc =read_inputfile('/user/home/lw13938/work/ciceroscm/inputfiles/AGAGE_conc.txt')
emi_input =read_inputfile('/user/home/lw13938/work/ciceroscm/inputfiles/AGAGE_em.txt')


def run_ciceroscm(gaspam_in, emi_input_in, df_ssp2_conc_in, conc_run=False, o3pert=0, taupert=0., Faci=0.):
    # NBVAL_IGNORE_OUTPUT
    scen = 'test'
    cscm_dir=CICEROSCM({
                "gaspam_data": gaspam_in,
                "emstart": 1750,  
                "conc_run":conc_run,
                "nystart": 1750,
                "nyend": 2024,
                "concentrations_data": df_ssp2_conc_in,
                "emissions_data": emi_input_in,
                "nat_ch4_data": df_nat_ch4,
                "nat_n2o_data": df_nat_n2o,
                "idtm":24,
            },o3pert=o3pert, taupert=taupert, Faci=Faci)

    # NBVAL_IGNORE_OUTPUT
    cscm_dir._run({
                "results_as_dict":True
            })
    out_df = to_df(cscm_dir)
    return out_df 

###############################################################################
# Run the ensemble
###############################################################################

forcing_keys = ["CH4", "STRAT_O3", "STRAT_H2O", "SO4_IND", "N2O", "TROP_O3", "SO4_DIR", "OC", "Total_forcing"]

df_ssp2_conc0 = df_ssp2_conc.copy()
emi_input0 = emi_input.copy()
for sp in csmsp:
    sp_in = check_hfc(sp)
    if sp == "CH3Br":
        df_ssp2_conc0[sp_in].values[:] = df_ssp2_conc[sp_in].values[0]
        emi_input0[sp_in].values[:] = emi_input[sp_in].values[0]
    else:
        df_ssp2_conc0[sp_in].values[:] = 0.
        emi_input0[sp_in].values[:] = 0.
gaspam_in = gaspam.copy()

# Run without emissions
df_out0  = run_ciceroscm(gaspam_in, emi_input0, df_ssp2_conc0, 
                                        conc_run=False)

years = df_out0['concentrations'].Year.values
years = years[np.isfinite(years)]
forcing_dict = {}
for key in forcing_keys:
    # Add each key to the dictionary
    forcing_dict[key] = np.zeros((len(years), nsamp))
forcing_dict["direct_forcing"] = np.zeros((len(years), nsamp))
forcing_dict["breakdown_products"] = np.zeros((len(years), nsamp))

# Use a sobol sequence to sample the parameter space
nvars = len(csmsp)*3 + 3
sobol = qmc.Sobol(d=nvars, scramble=True)
if int(fileindicator) == 1:
    seq = np.clip(sobol.random(nsamp), 1e-10, 1 - 1e-10)
else:
    seq = np.clip(sobol.fast_forward(int(fileindicator)-1).random(nsamp), 1e-10, 1 - 1e-10)

# Run the ensemble
for i in range(nsamp):
    gaspam_in = gaspam.copy()
    
    for sj,sp in enumerate(csmsp):
        sp_in = check_hfc(sp)
        # Hodnebrog et al. (2020) 90% uncertainty of 14% if TAU1 > 5, else 24%.
        # This equates to a sd of 8.51% and 14.59% respectively.
        # Also add the tropospheric adjustment to ALPHA so that and call all of this "direct RF".
        # That way I can separate out the breakdown products. 
        # CFC-11 and CFC-12 have their SARF to ERF conversions in Hodnebrog (2020). Others don't, use blanket number.
        if gaspam_in.loc[sp_in,"TAU1"] < 5:
            gaspam_in.loc[sp_in,"ALPHA"] = gaspam.loc[sp_in,"ALPHA"] * (1. + 0.1459* norm.ppf(seq[i,3*sj])) * \
                (1. + 0.079* norm.ppf(seq[i,3*sj+1]))
        else:
            if sp == "CFC-11":
                gaspam_in.loc[sp_in,"ALPHA"] = gaspam.loc[sp_in,"ALPHA"] * (1. + 0.0851* norm.ppf(seq[i,3*sj])) * \
                    (1.13 + 0.06* norm.ppf(seq[i,3*sj+1]))
            elif sp == "CFC-12":
                gaspam_in.loc[sp_in,"ALPHA"] = gaspam.loc[sp_in,"ALPHA"] * (1. + 0.0851* norm.ppf(seq[i,3*sj])) * \
                    (1.12 + 0.085* norm.ppf(seq[i,3*sj+1]))
            else:  
                gaspam_in.loc[sp_in,"ALPHA"] = gaspam.loc[sp_in,"ALPHA"] * (1. + 0.0851* norm.ppf(seq[i,3*sj])) * \
                    (1. + 0.079* norm.ppf(seq[i,3*sj+1]))
        
        # Some gases also have ERFs from breakdown products. Add these at this stage.
        if sp == "CFC-11":
            gaspam_in.loc[sp_in,"SARF_TO_ERF"] = 1.026 + 0.0032 * norm.ppf(seq[i,3*sj+2])
        elif sp == "CFC-12":
            gaspam_in.loc[sp_in,"SARF_TO_ERF"] = 1.015 + 0.0015 * norm.ppf(seq[i,3*sj+2])
        elif sp == "CFC-113":
            gaspam_in.loc[sp_in,"SARF_TO_ERF"] = 1.023 + 0.0022 * norm.ppf(seq[i,3*sj+2])
        elif sp == "HCFC-22":
            gaspam_in.loc[sp_in,"SARF_TO_ERF"] = 1.005 + 0.0005 * norm.ppf(seq[i,3*sj+2])
        elif sp == "CCl4":
            gaspam_in.loc[sp_in,"SARF_TO_ERF"] = 1.153 + 0.0174 * norm.ppf(seq[i,3*sj+2])
        else:
            gaspam_in.loc[sp_in,"SARF_TO_ERF"] = 1. # Other breakdown products are not considered
    
        
    # Functional fit for EESC to strat O3 SARF, CH4 lifetime and aerosols
    perto3 = -14.273 + 4.841 * norm.ppf(seq[i,-3]) #np.random.normal(loc=0, scale=2718.)
    taupert = 5.291e-5 + 2.547e-5 * norm.ppf(seq[i,-2])#np.random.normal(loc=0, scale=0.0014)
    Faci = 2.148e-6 + 4.436e-5 * norm.ppf(seq[i,-1])
    
    
    # Run without emissions (still need to run because of natural CH3Br)
    df_out0  = run_ciceroscm(gaspam_in, emi_input0, df_ssp2_conc0, 
                                            conc_run=False, o3pert=perto3, taupert=taupert, Faci=Faci)
           
    # Run with emissions        
    df_out = run_ciceroscm(gaspam_in, emi_input, df_ssp2_conc, 
                                        conc_run=False, o3pert=perto3, taupert=taupert, Faci=Faci)
    
    # Extract the forcing
    for key in forcing_keys:
        forcing_dict[key][:,i] = df_out["forcing"][key].dropna().values - df_out0["forcing"][key].dropna().values
    
    # Extract the direct forcing and breakdown products
    for sp in csmsp:
        sp_in = check_hfc(sp)
        dir_forcing_sp = df_out["concentrations"][sp_in].dropna().values * gaspam_in.loc[sp_in,"ALPHA"]
        forcing_dict["direct_forcing"][:,i] += dir_forcing_sp
        forcing_dict["breakdown_products"][:,i] += df_out["forcing"][sp_in].dropna().values - dir_forcing_sp

    if (i % 1000==0):
        print(f'Completed {i} members')

# Save to file
ds = xr.Dataset(
    {key: (["year", "sample"], data) for key, data in forcing_dict.items()},
    coords={"year": years, "sample": np.arange(nsamp)},
)

# Add global attributes
ds.attrs["title"] = f"Forcing Data of {gcase}"
ds.attrs["description"] = f"This dataset contains ERFs derived for {gcase} using an adaptation of the CICERO-SCM."
ds.attrs["history"] = f"Created on {pd.to_datetime('today').strftime('%Y/%m/%d')}."
ds.attrs["creator"] = "Luke Western"
ds.attrs["Conventions"] = "CF-1.8"

# Add variable-specific attributes
for key in forcing_dict.keys():
    ds[key].attrs["units"] = "W/m²"
    ds[key].attrs["long_name"] = f"{key.replace('_', ' ').title()} Effective Radiative Forcing"
    ds[key].attrs["description"] = f"Radiative forcing due to {key.replace('_', ' ').lower()}."

# Add coordinate-specific attributes
ds["year"].attrs["long_name"] = "Year"
ds["year"].attrs["standard_name"] = "time"
ds["sample"].attrs["long_name"] = "Sample Index"
ds["sample"].attrs["description"] = "Index of the Monte Carlo sample."

ds.to_netcdf(f"/user/home/lw13938/work/ciceroscm/RF_sensitivity/ensemble/outputs/ensemble_{gcase}_{fileindicator}.nc")
