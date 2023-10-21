import numpy as np
import pandas as pd
import fnmatch # for wildcard string searches


df_in = pd.read_csv("./icemodel.dat", sep = ' ' , names=(range(7)))
df_in.columns = ["z","scat_coef","abs_coef", "delta", "a1", "a2", "a3"]
new_scat_coef = np.ones_like(df_in["z"])*1E-10
new_abs_coef = np.ones_like(df_in["z"])*1E-10

dict = {"z" : df_in["z"], "scat_coef" : new_scat_coef, "abs_coef" : new_abs_coef, "delta" : df_in["delta"], "a1" : df_in["a1"], "a2" : df_in["a2"], "a3" : df_in["a3"]}
df_out  = pd.DataFrame(dict)

df_out.to_csv("icemodel_new.dat", header=False, index=False, sep = " ", float_format="%.10f")
