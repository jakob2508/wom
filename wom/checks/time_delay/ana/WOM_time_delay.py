import numpy as np
from scipy.interpolate import PchipInterpolator

def convolution(x, func1, func2):
    # discrete numerical convolution
    index = np.arange(len(x))
    conv = sum([func2[(index-i)]*func1[i] for i in index])
    return conv

def exponential(x, arg):
    # exponential function
    tau = arg
    f = np.exp(-x/tau)
    return f/(np.sum(f)*(x[1]-x[0]))

def moffat_king(x, par):
    # moffat-king function
    sigma, gamma = par
    y = x/(2*np.pi*sigma**2) * (1 - 1/gamma) * (1 + x**2/(2*gamma * sigma**2))**(-gamma)
    return y

def moffat_king_exponential(x, par, arg):
    # convolution of moffat-king and exponential
    A, sigma, gamma = par
    tau = arg
    func1 = moffat_king(x, (sigma, gamma))
    func2 = exponential(x, tau)
    y = convolution(x, func1, func2)
    y = y/(np.sum(y) * (x[1]-x[0]))
    y *= A
    return y

def get_time_sort_zpos(hit_time, hit_zpos, hit_pmt):
    # Bins arrival times of hits given their z hit position
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
    for zz in np.arange(wom_length/2-wzs, -wom_length/2, -wzs): # binning is from +-37 cm not +-38 cm
        zz = np.round(zz, 2) # zz rounded to 2 digits precision

        if zz == wom_length/2-wzs: # edge case hit at +38cm assigned to distribiution at +37cm
            zpos_mask = np.logical_and((hit_zpos>(zz-wzs/2)),(hit_zpos<=zz+wzs)).astype('bool')
        elif zz == -wom_length/2+wzs: # edge case hit at -38cm assigned to distribiution at -37cm
            zpos_mask = np.logical_and((hit_zpos>(zz-wzs)),(hit_zpos<=(zz+wzs/2))).astype('bool')
        else:
            zpos_mask = np.logical_and((hit_zpos>(zz-wzs/2)),(hit_zpos<=(zz+wzs/2))).astype('bool')
        time_sort.append(hit_time[zpos_mask])
        n_samples.append(np.sum(zpos_mask))
        new_order.append(order[zpos_mask])
    return time_sort, n_samples, new_order

def get_new_time(time, time_delay, new_order):
    # Combines arrival time and time delay and puts them back into the original order before the histogramming of data
    mask_sort = np.argsort(np.concatenate(new_order).flatten())
    new_time = []
    for j in range(len(time)):
        new_time.append(time[j]+time_delay[j])

    new_time = np.concatenate(new_time).flatten()
    sort_time = new_time[mask_sort]
    return sort_time

def get_pdf(moffat_king_exponential, x, td_zpos, td_params, tau):
    # Returns array of discretized pdfs following a convolution between a Moffat-King and an exponential distribution
    pdf = []
    bin_width = 0.02
    for zpos, tdp in zip(td_zpos,td_params):
        # calculate minimal geomtric time offset to be added to moffat-king time delay
        t_offset = (wom_length/2-zpos)/(c/n_ref)*1E9
        # corresponding bin offset
        bin_offset = round(t_offset/bin_width)
        y = moffat_king_exponential(x, par=tdp, arg=tau)
        # shift Moffat-King distribution
        y = np.roll(y, bin_offset)
        pdf.append(y)
    return np.array(pdf)

def get_num_cdf(pdf):
    # Numerical cdf given a pdf
    cdf = np.cumsum(pdf, axis = -1)
    cdf /= cdf[:,-1, np.newaxis]
    return cdf
  
def get_inv_cdf_interpol(x,pdf):
    # Uses spline interpolation to build the inverse cdf
    cdf = get_num_cdf(pdf)
    inv_cdf_interpol = []
    for c in cdf:
        inv_cdf_interpol.append(PchipInterpolator(c, x))
    return inv_cdf_interpol

def get_time_delay(x, pdf, n_samples):
    # Applies inverse cdf sampling to pull random time delay from pdf
    inv_cdf = get_inv_cdf_interpol(x, pdf)
    time_delay = []
    for j, nn in enumerate(n_samples):
        rand_cdf = np.random.uniform(0,1,size=nn)
        time_delay.append(inv_cdf[j](rand_cdf))
    return time_delay

def get_time_delay_loop(x, pdf, n_samples):
    # Applies inverse cdf sampling to pull random time delay from pdf, modified for 'loop' simulation mode
    inv_cdf = get_inv_cdf_interpol(x, pdf)
    rand_cdf = np.random.uniform(0,1,size=n_samples)
    time_delay = inv_cdf[0](rand_cdf)
    return time_delay

def sim_time_delay(hit_time, hit_zpos, hit_pmt, sim_mode = 'serial'):

    # Time binning [0, 20] ns in 1001 steps, Used to construct the pdf, cdf and for sampling.
    x = np.linspace(0,20, 1001, endpoint=True)

    # Load time delay pdf parametrization from file
    time_delay_dist_data = np.genfromtxt('td_140.0', delimiter='')

    td_zpos = time_delay_dist_data.T[0]
    td_params = time_delay_dist_data.T[1:].T

    # Get array of numerical pdf values of the z dependent time delay distribution
    pdfs = get_pdf(moffat_king_exponential, x, td_zpos, td_params, tau)


    # In the 'serial' simulation mode, all hits get histogrammed and the time delays are pulled from the pdfs based on the histogram entries.
    # In the 'loop' simulation mode every hit pulls exactly 1 time delay from the respective pdf distribution of the z hit position.
    # For a file of about 5500 hits, the 'serial' mode takes 0.01 s while the 'loop' mode takes 0.78 s, alomost 80 times longer. It is 
    # thus recommended to use the 'serial' mode. The results between the different modes should vanish for large statistical sample as a result
    # of the law of large numbers.
    if sim_mode == 'serial':
        time_sort, n_samples, new_order = get_time_sort_zpos(hit_time, hit_zpos, hit_pmt)
        time_delay = get_time_delay(x, pdfs, n_samples)
        time_new = get_new_time(time_sort, time_delay, new_order)
        return time_new
    elif sim_mode == 'loop':
        # Prepare the z hit position data: Flipping PMT position, rounding and shifting edge hits
        hit_zpos = np.where(hit_pmt == 0, -hit_zpos, hit_zpos)
        hit_zpos = np.round(hit_zpos, 2)
        # The simulation of time delays was done between [-37,37] cm in z while the tube is simulated from [-38,38] cm
        # For simplicity we treat hits in +-38 cm as if they hit +-37 cm which does not make a big difference.
        hit_zpos = np.where(hit_zpos == 0.38, 0.37, hit_zpos)
        hit_zpos = np.where(hit_zpos == -0.38, -0.37, hit_zpos)

        time_delay = []
        for i in range(len(hit_time)):
            hz = hit_zpos[i]
            td = get_time_delay_loop(x, pdfs[td_zpos == hz], 1)
            time_delay.append(td)
        time_delay = np.array(time_delay).flatten()
        time_new = hit_time + time_delay
        return time_new

# Some constants
global c
global n_ref
global tau
global wom_length
global wom_z_sampling

c = 299792458
n_ref = 1.46
tau = 1.6
wom_length = 0.76
wom_z_sampling = 0.01
