#Plot radiative forcing all components in the CICERO SCM
#Not all components have radiative forcing output.

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 12

#Read results:
outdir = '/div/no-backup/users/ragnhibs/ciceroscm/scripts/output_test_forcing/'
scenlist = ['test_v_sun','test_u_sun']


fig, axs = plt.subplots(nrows=1, ncols=1,sharex=True,figsize=(12,8))

fig.suptitle('CICERO SCM simulation')
for scen in scenlist:
    print(scen)
    
    df_rf=pd.read_csv(outdir+'/'+ scen +'_forc.txt', sep='\t', index_col=0)
    print(df_rf)

    #Plot first 16 components:
    complist = ['Total_forcing']
    
    for i,c in enumerate(complist):
        print(i)
        print(c)
        print(df_rf.columns)
        comp = c 
        print(comp)
        df_rf[comp].plot(ylabel='RF [Wm$^{-2}$ ]',ax=axs,label=scen)
        axs.set_title(comp)
        axs.legend()
        axs.axhline(y=0,color='k',linestyle=':',linewidth=0.5)
    
 

plt.show()
exit()
