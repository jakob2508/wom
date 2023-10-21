import numpy as np
import pandas as pd
import fnmatch # for wildcard string searches
import matplotlib.pyplot as plt

file = pd.read_csv("../sim/logs/log.txt", sep = ' ' , names=(range(13))) # load table from log
file = file.drop(file[file[0]!="HIT"].index) # drops lines without HIT
file = file.dropna(how='all', axis = 1) # drops lines with NaN
file.columns = ["HIT","str","om_pmt", "time", "wavelength", "pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z"] # column names
wvl = file["wavelength"].values
wvl = np.array(wvl, dtype = float)

# plot histogram of wavelength
fig, ax = plt.subplots(1,1)
ax.hist(wvl, bins = 50, histtype="step")
ax.set_xlabel("wavelength [nm]")
ax.set_ylabel("entries")
plt.tight_layout()
plt.savefig('./plots/check_eff_area.png')
