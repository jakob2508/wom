import numpy as np
import matplotlib.pyplot as plt
from WOM_time_delay import *
import pandas as pd

debug1 = 0 # plots pdf and cdf and 1000 randomly drawn samples
debug2 = 1 # plots histogram of PMT up/down and the combined distribution
debug3 = 1 # histograms arrival time of photons and time delayed photon distribution

str_flash = 1 # str number flasher
om_flash = 1 # om number flasher
str_wom = 2 # str number wom
om_wom = 1 # om number wom

if debug1:
    # time binning of pdf, [0, 20] ns
    x = np.linspace(0,20, 1001, endpoint=True)

    # import time delay parameters
    time_delay_dist_data = np.genfromtxt('td_140.0', delimiter='')

    td_zpos = time_delay_dist_data.T[0]
    td_params = time_delay_dist_data.T[1:].T

    # uniform binning for every z position
    n_samples = np.array(np.ones(len(td_zpos))*1000, dtype=int)
    pdfs = get_pdf(moffat_king_exponential, x, td_zpos, td_params, tau)
    cdfs = get_num_cdf(pdfs)
    inv_cdfs = get_inv_cdf_interpol(x, pdfs)
    inv_sampling = get_time_delay(x, pdfs, n_samples)

    for i in range(len(td_params)):
        fig, ax = plt.subplots(1,1)
        ax2 = ax.twinx()
        ax.plot(x, pdfs[i], color = 'C0')
        ax2.plot(x, cdfs[i], color = 'C1')
        ax.hist(inv_sampling[i], bins = 50, density=True, range = (0,20), histtype = 'step', color = 'grey')
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('pdf')
        ax2.set_ylabel('cdf')
        ax2.set_ylim((0,1))
        plt.savefig('./plots/test_pdf_cdf_sample_z={:.0f}cm.png'.format(td_zpos[i]*100))
        plt.close()

# load and prepare data
file = pd.read_csv("./../ang_sens/sim/logs/log_theta=90.txt", sep = ' ' , names=(range(13))) # load table from log
file = file.drop(file[file[0]!="HIT"].index) # drops lines without HIT
file = file.dropna(how='all', axis = 1) # drops lines with NaN
file.columns = ["HIT","str","om_pmt", "time", "wavelength", "pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z"] # column names

hit_pmt = file["om_pmt"].str.contains("_1").values.astype(int)
hit_zpos = file["pos_z"].values.astype(float)
hit_time = file["time"].values.astype(float)

if debug2:
    time_sort, n_samples, new_order = get_time_sort_zpos(hit_time, hit_zpos, hit_pmt)
    hist_y, bins = np.histogram(hit_zpos, bins = 77, range=(-0.38,0.38), density=True)
    hist_x = (bins[1:]+bins[:-1])/2 * 100

    hit_zpos_pmt0 = hit_zpos[hit_pmt == 0]
    hit_zpos_pmt1 = hit_zpos[hit_pmt == 1]
    hit_zpos_pos = np.concatenate([-hit_zpos_pmt0, hit_zpos_pmt1]) # PMT 0/DOWN -> flip

    fig, ax = plt.subplots(1,2, figsize = (10,4))
    ax[0].hist(hit_zpos*100, bins = 77, color = 'grey', range = (-38.5, 38.5), label = r'$H(z)$')
    ax[0].hist(hit_zpos_pmt0*100, bins = 77, histtype='step', lw = 3, range = (-38.5, 38.5), label = r'$H_{down}(z)$')
    ax[0].hist(hit_zpos_pmt1*100, bins = 77, histtype='step', lw = 3, range = (-38.5, 38.5), label = r'$H_{up}(z)$')

    ax[1].hist(hit_zpos_pos*100, bins = 75, histtype='step', lw = 3, color = 'grey', range = (-37.5, 37.5), label = r'$H_{down}(z) + H_{up}(-z), rebin$')
    ax[1].step(np.arange(37, -38, -1), n_samples, lw = 3, color = 'red', ls = '-.', label = r'resort algorithm', where = 'mid')
    
    for i in range(2):
        ax[i].set_xlabel('z position [cm]')
        ax[i].set_ylabel('entries')
        ax[i].legend()
    plt.tight_layout()
    plt.savefig('./plots/time_delay_zpos_split.png')

# simulate time delay
time_new_serial = sim_time_delay(hit_time, hit_zpos, hit_pmt, sim_mode='serial')
time_new_loop = sim_time_delay(hit_time, hit_zpos, hit_pmt, sim_mode='loop')

if debug3:
    fig, ax = plt.subplots(1,1)
    ax.hist(hit_time, bins = 100, range = (14,35), histtype='step', density=True, label=r'$t_{in}$')
    ax.hist(time_new_serial, bins = 100, range = (14,35), histtype='step', density=True, label=r'$t_{out}^{serial}$')
    ax.hist(time_new_loop, bins = 100, range = (14,35), histtype='step', ls = '--', density=True, label=r'$t_{out}^{loop}$')
    ax.set_xlabel('Time [ns]')
    ax.set_ylabel('Normalized Entries')
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()
    plt.savefig('./plots/time_delay_final.png')