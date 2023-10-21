import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad
from scipy.special import expi

def convolution(x, func1, func2, D):
    # discrete numerical convolution
    y1 = func1(x)
    y2 = func2(x, D)
    index = np.arange(len(y1))
    conv = sum([y1[(index-i)]*y2[i] for i in index])
    
    return conv

def exponential(t):
    # exponential time delay of paint
    f = 1/tau * np.exp(-t/tau)
    return f

def geometry(t, D):
    # geometric time delay
    xx = t * c/n_ref
    f = D**2/xx**2 * np.exp(-xx/L_att) * np.heaviside(xx-D, 1)
    norm = (n_ref * D)/c * (np.exp(-D/L_att) + D/L_att * expi(-D/L_att))
    f = f/norm
    f = np.where(np.isnan(f) == True, 0, f)
    return f

def combined(t, tp, D):
    return exponential(t-tp)*geometry(t, D)

def conv(t_int, D):
    print(t_int)
    #i = []
    #for tp in t_int:
    c, foobar = quad(combined, 0, 1E-6, args = (t_int, D)) 
    #    i.append(c)
    print(c)
    return np.array(c)

def get_time_sort_zpos(hit_time, hit_zpos, hit_pmt):
    order = np.arange(len(hit_time))
    time_sort = [] # contains photon time sorted for z bin
    n_samples = [] # number of hits in that z bin
    new_order = [] # save the order in which we re-shuffle the time array
    wzs = wom_z_sampling
    # those hits that are assinged to PMT 0 (DOWN) are treated as if they hit on PMT 1 but with the sign of the z position flipped
    # e.g. a photon with the signature (t, z=-37 cm, PMT = 0) == (t, z=37 cm, PMT = 1) because used PMT 1 (UP) for the time delay
    # characterisation
    hit_zpos = np.where(hit_pmt == 0, -hit_zpos, hit_zpos)
    # group/histogram hits in bins of 1 cm with special care for the edge cases
    for zz in np.arange(wom_length/2, -wom_length-wzs, -wzs):
        if zz == wom_length/2: # edge case hit at +38cm assigned to distribiution at +37cm
            zpos_mask = np.logical_and((hit_zpos>(zz-wzs/2)),(hit_zpos<(zz+wzs))).astype('bool')
        elif zz == -wom_length/2: # edge case hit at -38cm assigned to distribiution at -37cm
            zpos_mask = np.logical_and((hit_zpos>(zz-wzs)),(hit_zpos<(zz+wzs/2))).astype('bool')
        else:
            zpos_mask = np.logical_and((hit_zpos>(zz-wzs/2)),(hit_zpos<(zz+wzs/2))).astype('bool')
        time_sort.append(hit_time[zpos_mask])
        n_samples.append(np.sum(zpos_mask))
        new_order.append(order[zpos_mask])
    return time_sort, n_samples, new_order

def get_new_time(time, time_delay, n_samples, new_order):
    mask_sort = np.argsort(new_order)
    new_time = []
    for j in range(len(n_samples)):
        new_time.append(time[j]+time_delay[j])

    new_time = np.array(new_time, dtype='object')
    sort_time = new_time[mask_sort]
    return sort_time

def get_pdf(moffat_king_exponential, x, td_zpos, td_params, tau):
    pdf = []
    bin_width = 0.02
    for zpos, tdp in zip(td_zpos,td_params):
        # calculate minimal geomtric time offset to be added to moffat-king time delay
        t_offset = zpos/(c/n_ref)*1E9
        # corresponding bin offset
        bin_offset = round(t_offset/bin_width)
        y = moffat_king_exponential(x, par=tdp, arg=tau)
        # shift Moffat-King distribution
        y = np.roll(y, bin_offset)
        pdf.append(y)
    return np.array(pdf)

def get_num_cdf(pdf):
    cdf = np.cumsum(pdf, axis = -1)
    cdf /= cdf[:,-1, np.newaxis]
    return cdf
  
def get_inv_cdf_interpol(x,pdf):
    cdf = get_num_cdf(pdf)
    inv_cdf_interpol = []
    for c in cdf:
        inv_cdf_interpol.append(PchipInterpolator(c, x))
    return inv_cdf_interpol

def get_time_delay(x, pdf, n_samples):
    # applies inverse cdf sampling to pull random time delay from pdf
    inv_cdf = get_inv_cdf_interpol(x, pdf)
    time_delay = []
    for j, nn in enumerate(n_samples):
        rand_cdf = np.random.uniform(0,1,size=nn)
        time_delay.append(inv_cdf[j](rand_cdf))
    return time_delay

global c
global n_ref
global tau
global wom_length
global wom_z_sampling
global L_att

c = 299792458
n_ref = 1.46
tau = 1.6E-9
wom_length = 0.76
wom_z_sampling = 0.01
L_att = 3


# strategy
# sort times according to pmt and z pos
# initialise moffat-king-exponential distribution once
# interpolate distribution
# sample via inverse transform sampling
# numerically calculate inverse function by looking for x for which cdf(x) = p, loop over p
# consider using splines instead of custom minimzation
# don't forget tome shift!!!!
# PMT selection !!!
# change from 0, 76 cm to -33 to 33 cm
# filter for PMT number first
# add lower PMT (PMT = 0 ?) with inverted zpos to arrays
# get time delay for ordered z range
# re-sort array according to PMT ID and occurance in log file

# test strategy:
# test sorting input hits by injecting mono bins, sum of all hist entries = len(input)
# test moffat_king pdf and cdf visually - DONE!
# test inverse sampling with gaussian and compare to scipy.special - DONE!