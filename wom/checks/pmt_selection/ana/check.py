import numpy as np
import pandas as pd
import fnmatch # for wildcard string searches
import matplotlib.pyplot as plt

from plthelper import *

str_flash = 1 # str number flasher
om_flash = 1 # om number flasher
str_wom = 2 # str number wom
om_wom = 1 # om number wom

# postion of WOM (1,1) in the IceCube coordinate frame
# nominal position
wom_x0 = 1
wom_y0 = 0
wom_z0 = 1948.07-2300
wom_pos0 = np.array([wom_x0,wom_y0,wom_z0])

# actual position
wom_x = 1
wom_y = 0
wom_z = 1948.07-2300
wom_pos = np.array([wom_x,wom_y,wom_z])

wom_dpos = wom_pos-wom_pos0

z_range = np.array([-.30,-.25,-.20,-.15,-.10,-.05,.00,.05,.10,.15,.20,.25,.30])

NUM_PH = []
NUM_HIT = []
NUM_PMT1 = []
FLA_dpos = []
FLA_theta = []
FLA_phi = []
FLA_opening_angle = []
FLA_length = []
HIT = []

for z in z_range:

    file = pd.read_csv("../sim/logs/log_z={:.2f}.txt".format(z), sep = ' ' , names=(range(13))) # load table from log
    num = file[file[0].str.contains("photons").astype("boolean")]
    num = num.dropna(how='all', axis = 1) # drops lines with NaN
    num.columns = ["ph_txt","ph_num","hit_txt", "hit_num"]
    flash1 = file[file[0].str.contains("Flasher").astype("boolean")]
    flash2 = file[file[1].str.contains("flasher").astype("boolean")]
    file = file.drop(file[file[0]!="HIT"].index) # drops lines without HIT
    file = file.dropna(how='all', axis = 1) # drops lines with NaN
    file.columns = ["HIT","str","om_pmt", "time", "wavelength", "pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z"] # column names
    dirz = file["pos_y"].values
    dirz = np.array(dirz, dtype = float)

    str_mask = file["str"].str.contains(str(str_wom)).values # string mask
    om_mask = file["om_pmt"].str.contains(str(om_wom)).values # OM mask
    pmt_mask = file["om_pmt"].str.contains("_1").values
    mask = (str_mask*om_mask).astype(bool)
    mask2 = (mask*pmt_mask).astype(bool)

    hit_x = file["pos_x"].values[mask].astype(float)
    hit_y = file["pos_y"].values[mask].astype(float)
    hit_z = file["pos_z"].values[mask].astype(float)
    hit_pos = np.array([hit_x,hit_y,hit_z])

    fla_pos_x = float(flash1.iloc[0][3].strip(","))
    fla_pos_y = float(flash1.iloc[0][4].strip(","))
    fla_pos_z = float(flash1.iloc[0][5].strip(","))
    fla_phi = float(flash1.iloc[1][11])
    fla_opening_angle = float(flash2.iloc[1][5])/2 * np.pi/180
    fla_theta = float(flash2.iloc[2][5])
    fla_theta_cor = -0.2
    fla_theta -= fla_theta_cor
    fla_pos = np.array([fla_pos_x, fla_pos_y, fla_pos_z])
    fla_dpos = fla_pos-wom_pos0
    #fla_length = np.sqrt(np.sum((wom_origin-fla_origin)**2))
    fla_length = np.sqrt(np.sum((np.mean(hit_pos,axis=1)-fla_dpos)**2))

    NUM_PH.append(np.sum(num["ph_num"].values.astype(int)))
    NUM_HIT.append(mask.sum())
    NUM_PMT1.append(mask2.sum())
    HIT.append(hit_pos)
    FLA_dpos.append(fla_dpos)
    FLA_theta.append(fla_theta)
    FLA_phi.append(fla_phi)
    FLA_opening_angle.append(fla_opening_angle)
    FLA_length.append(fla_length)

NUM_PH = np.array(NUM_PH)
NUM_HIT = np.array(NUM_HIT)
NUM_PMT1 = np.array(NUM_PMT1)



rho = 0.057 # 110 mm cylinder outer diameter # 0.057 in om.conf
height = 2*0.38 # 640 mm inner cylinder height # 0.38 in om.conf
color = ["C0", "C1", "C2", "C3", "C4"]
skip = 10
wd = 1
j = 0

if 0:
    fig, ax = plt.subplots(1, 1, figsize = (4,4), subplot_kw={'projection': '3d'})
    ax.set_box_aspect((1,1,1))
    ax = plot_3D_cylinder(ax, rho, height, x0 = wom_dpos[0], y0 = wom_dpos[1], z0 = wom_dpos[2], color = "k")
    for i in range(0,len(z_range),3):
        ax = plot_3D_sphere(ax, 0.01, x0 = FLA_dpos[i][0], y0 = FLA_dpos[i][1], z0 = FLA_dpos[i][2], color = "red")
        ax = plot_3D_cone(ax, FLA_dpos[i], phi = FLA_phi[i], theta = FLA_theta[i], length = FLA_length[i], opening_angle= FLA_opening_angle[i], color = color[j])
        ax.scatter(HIT[i][0][::skip], HIT[i][1][::skip], HIT[i][2][::skip], alpha = 1, color = color[j], s = 1)
        j += 1


    # ax.text(wom_x+rho,wom_y,wom_z, "[0,0,0]", zdir=(0,1,0))
    # ax.text(wom_x,wom_y+rho,wom_z, "[0,90,0]", zdir=(1,0,0))
    # ax.text(wom_x-rho,wom_y,wom_z, "[0,180,0]", zdir=(0,-1,0))
    # ax.text(wom_x,wom_y-rho,wom_z, "[0,270,0]", zdir=(-1,0,0))
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.set_xlim(wom_dpos[0]-wd,wom_dpos[0]+wd)
    ax.set_ylim(wom_dpos[1]-wd,wom_dpos[1]+wd)
    ax.set_zlim(wom_dpos[2]-wd,wom_dpos[2]+wd)
    ax.set_xticks([wom_dpos[0]-wd,wom_dpos[0],wom_dpos[0]+wd])
    ax.set_yticks([wom_dpos[1]-wd,wom_dpos[1],wom_dpos[1]+wd])
    ax.set_zticks([wom_dpos[2]-wd,wom_dpos[2],wom_dpos[2]+wd])
    plt.tight_layout()
    plt.savefig("./plots/sim_setup.png")

if 0:
    fig, ax = plt.subplots(1, 1)
    ax.hist(HIT[0][2], range=[-0.32,0.32], bins = 20, histtype="step", color = "C0", label = r"-45$\deg$")
    ax.hist(HIT[-1][2], range=[-0.32,0.32], bins = 20, histtype="step", color = "C1", label = r"+45$\deg$")
    ax.set_xlabel('z hit position')
    ax.set_ylabel('entries')
    ax.legend()
    plt.tight_layout()
    plt.savefig("./plots/hist_posz_cones.png")
    plt.show()

if 1:
    z_range = (height/2 - z_range)*100
    NUM_PMT0 = NUM_HIT-NUM_PMT1
    R_PMT0 = NUM_PMT0/NUM_HIT * 100
    R_PMT1 = NUM_PMT1/NUM_HIT * 100
    UR_PMT0 = np.sqrt(NUM_PMT0)/NUM_HIT * 100
    UR_PMT1 = np.sqrt(NUM_PMT1)/NUM_HIT * 100
    fig, ax = plt.subplots(1, 1)
    ax.errorbar(x=z_range,y=R_PMT0, yerr=UR_PMT0, label = r'PMT$_{down}$', color = 'k', capsize=3)
    ax.errorbar(x=z_range,y=R_PMT1, yerr=UR_PMT1, label = r'PMT$_{up}$', color = 'grey', capsize=3)
    ax.axhline(100, color = 'red')
    ax.set_xlabel(r'z position [cm]')
    ax.set_ylabel('normalized counts')
    ax.legend()
    ax.set_xlim(0, height*100)
    ax.set_ylim(0,110)
    plt.tight_layout()
    plt.savefig("./plots/check_pmt_selection.png")
    plt.show()