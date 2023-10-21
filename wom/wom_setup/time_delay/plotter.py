import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator
import numpy as np
import pickle
from iminuit import Minuit

generate = 1
debug1 = 0 # plots the tts for different z positions
debug2 = 0 # plots the combined tts for all z positions
debug3 = 0 # plots the quantiles, max, median, FWHM as a function of the z position
debug4 = 0 # plot fit parameter for Moffat-King exponential fit as function of distance
debug5 = 0 # plot pdf as function of distance between 0 and 74 cm
debug6 = 0

# returns median, maximum, FWHM, quantiles of histogramed data
def hist_charac(hist_x,hist_y):
    dx = hist_x[1] - hist_x[0]
    hist_max = np.max(hist_y)
    ind_max = np.argmax(hist_y)
    t_max = hist_x[ind_max]
    t_lhm = float(hist_x[np.argwhere(hist_y >= hist_max/2)[0]])
    t_uhm = float(hist_x[np.argwhere(hist_y >= hist_max/2)[-1]])
    t_50 = hist_x[np.argsort(np.abs(np.cumsum(hist_y)*dx-0.5))[0]]
    t_lq90 = hist_x[np.argsort(np.abs(np.cumsum(hist_y)*dx-0.1))[0]]
    t_uq90 = hist_x[np.argsort(np.abs(np.cumsum(hist_y)*dx-0.9))[0]]
    t_lq95 = hist_x[np.argsort(np.abs(np.cumsum(hist_y)*dx-0.05))[0]]
    t_uq95 = hist_x[np.argsort(np.abs(np.cumsum(hist_y)*dx-0.95))[0]]    
    return t_max, t_50, t_lq90, t_uq90, t_lq95, t_uq95, t_lhm, t_uhm

def convolution(x, hist_y, fold):
    # convolutes hist_y and fold 
    conv = sum([fold[(x-i)]*hist_y[i] for i in range(len(hist_y))])
    return conv

def exponential(x, arg):
    # exponential folding function
    tau = arg
    f = np.exp(-x/tau)
    return f/(np.sum(f)*(x[1]-x[0]))

def hyperbolic(x, arg):
    z = arg
    t_offset = z*1E-2/(c/n_ref)*1E9
    y = np.where(x>=t_offset, 1/x**2, 0)
    return y

def moffat_king(x, par):
    sigma, gamma = par
    y = x/(2*np.pi*sigma**2) * (1 - 1/gamma) * (1 + x**2/(2*gamma * sigma**2))**(-gamma)
    return y

def hyperbolic_exponential(x, par, arg):
    A = par
    x, z, tau = arg
    t_offset = z*1E-2/(c/n_ref)*1E9
    func1 = hyperbolic(x, z)
    func2 = exponential(x, tau)
    y = convolution(np.arange(len(x)), func1, func2)
    y = np.where(x>=t_offset, y, 0)
    y = y/(np.sum(y) * (x[1]-x[0]))
    y *= A
    return y

def moffat_king_exponential(x, par, arg):
    A, sigma, gamma = par
    tau = arg
    #t_offset = z*1E-2/(c/n_ref)*1E9
    func1 = moffat_king(x, (sigma, gamma))
    func2 = exponential(x, tau)
    y = convolution(np.arange(len(x)), func1, func2)
    #y = np.where(x>=t_offset, y, 0)
    y = y/(np.sum(y) * (x[1]-x[0]))
    y *= A
    return y

def loss(par, args):
    # objective function for scipy.minimize
    bin_width = 0.1
    x, y_true, z, tau = args
    dof = len(x) - len(par)
    t_offset = z*1E-2/(c/n_ref)*1E9
    bin_offset = round(t_offset/bin_width) 
    y_true = np.roll(y_true, -bin_offset)
    #y_reco = moffat_king(x, par)
    #y_reco = hyperbolic(x, par, arg)
    #y_reco = hyperbolic_exponential(x, par, arg=(x,z,tau))
    y_reco = moffat_king_exponential(x, par, arg=tau)
    loss = np.sqrt(np.sum((y_reco-y_true)**2))/dof
    return loss

'''
def obj(par):
    # objective function for iminuit.Minuit
    bin_width = 0.1
    dof = len(x) - len(par)
    y_true = fold_y
    t_offset = z*1E-2/(c/n_ref)*1E9
    bin_offset = round(t_offset/bin_width) 
    y_true = np.roll(y_true, -bin_offset)
    y_reco = moffat_king_exponential(x, par, tau)
    loss = np.sqrt(np.sum((y_reco-y_true)**2))/dof
    return loss
'''

def moving_average(a, n=3, zero_padding = False, const_padding = False):
    # moving average filter used to smoothen the fit parameter dependency on the hit position z
    if zero_padding:
        ind = np.arange(n-1)
        a = np.insert(a, ind, np.zeros(n-1))
        a = np.roll(a, -int((n-1)/2))
    if const_padding:
        if n%2 != 1:
            ind1 = np.arange((n-1)/2).astype(int)
            ind2 = -np.arange(1,(n-1)/2).astype(int)
            a = np.insert(a, ind1, np.ones(int((n/2)))*a[0])
            a = np.insert(a, ind2, np.ones(int((n-1)/2))*a[-1])
        else:
            ind1 = np.arange((n-1)/2).astype(int)
            ind2 = -np.arange(1,(n+1)/2).astype(int)
            a = np.insert(a, ind1, np.ones(int((n-1)/2))*a[0])
            a = np.insert(a, ind2, np.ones(int((n-1)/2))*a[-1])
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
    

# loading data
file = open("./files/sim_delay_1cm_sampling.pkl","rb")
z_range, distances1, distances2 = pickle.load(file)
file.close()

# which pmt?
pmt = "up" #"up", "down"

distances1 = np.array(distances1, dtype = object)
distances2 = np.array(distances2, dtype = object)


# ALGO simulation saves traveled distance in [cm] and not time in [ns]
times1 = distances1 * 0.01 * 1.47/(299792458) * 1E9
times2 = distances2 * 0.01 * 1.47/(299792458) * 1E9

c = 299792458
n_ref = 1.46

wom_length = 0.76

FIT = [] # fitting parameters
T = []
Y = []
T_50 = []
T_max = []
T_lq90 = []
T_lq95 = []
T_uq90 = []
T_uq95 = []
T_lhm = []
T_uhm = []

if generate:

    # prepare data: loop over files and save median, maximum, quantiles, FWHM
    for il in range(len(z_range)):
        if pmt == "up":
            t = times1[il]
        elif pmt == "down":
            t = times2[il]
        # histogram timing distribution
        hist_y, bins = np.histogram(t, bins = 200, density=True, range = (0,20))
        hist_x = (bins[1:]+bins[:-1])/2

        # paint folding function
        tau = 1.6
        paint = exponential(hist_x, tau)
        paint /= np.sum(paint)

        # folding
        fold_y=convolution(np.arange(len(hist_x)), hist_y, paint)
        fold_y = fold_y/(np.sum(fold_y) * (hist_x[1] - hist_x[0]))

        t_max, t_50, t_lq90, t_uq90, t_lq95, t_uq95, t_lhm, t_uhm = hist_charac(hist_x, fold_y)
        
        #res = minimize(loss, x0 = ((0.5,0.5,0.5)))
        if pmt == 'up':
            res = minimize(loss, x0 = ((0.5,0.5,0.5)), args = [hist_x, fold_y, z_range[il], tau])
        elif pmt == 'down':
            res = minimize(loss, x0 = ((0.5,0.5,0.5)), args = [hist_x, fold_y, z_range[-1-il], tau])

        '''
        x, y_true, z, tau = hist_x, fold_y, z_range[il], tau
        ini_A, ini_sigma, ini_gamma = 0.5, 0.5, 0.5
        opt = Minuit(obj, (ini_A, ini_sigma, ini_gamma), name = ('A', 'sigma', 'gamma'))
        opt.tol = 1E-1
        opt.errors = (1E-2, 1E-2, 1E-2)
        opt.errordef = Minuit.LEAST_SQUARES
        opt.limits = [(0,None),(None,None),(0,None)]
        opt.strategy = 2
        
        opt.simplex()

        opt.fixed['A'] = False
        opt.fixed['sigma'] = True
        opt.fixed['gamma'] = False
        opt.migrad()

        opt.fixed['A'] = True
        opt.fixed['sigma'] = False
        opt.fixed['gamma'] = True
        opt.migrad()

        opt.fixed['A'] = False
        opt.fixed['sigma'] = False
        opt.fixed['gamma'] = False
        opt.migrad()

        FIT.append(opt.values[:])
        '''
        FIT.append(res.x)
        T.append(hist_x)
        Y.append(fold_y)
        T_50.append(t_50)
        T_max.append(t_max)
        T_lq90.append(t_lq90)
        T_uq90.append(t_uq90)
        T_lq95.append(t_lq95)
        T_uq95.append(t_uq95)
        T_lhm.append(t_lhm)
        T_uhm.append(t_uhm)

        if il == 18: # save the pdf of central illumination at 38 cm (il = 18)
            file = open("./pdf_wom_z=38cm.pkl","wb")
            pickle.dump([hist_x, fold_y], file)
            file.close()

    FIT = np.array(FIT)
    T = np.array(T)
    Y = np.array(Y)
    T_50 = np.array(T_50)
    T_max = np.array(T_max)
    T_lq90 = np.array(T_lq90)
    T_uq90 = np.array(T_uq90)
    T_lq95 = np.array(T_lq95)
    T_uq95= np.array(T_uq95)
    T_lhm = np.array(T_lhm)
    T_uhm= np.array(T_uhm)

    z_pos = wom_length/2 - z_range*1E-2

    filename = 'td_140.0_'+pmt+'.txt'
    data = np.vstack((z_pos, FIT.T[0], FIT.T[1], FIT.T[2])).T
    np.savetxt(filename, data, delimiter=' ', newline='\n', fmt = ('%.2f', '%.8f','%.8f','%.8f'))


FIT_up = np.genfromtxt('td_140.0_up.txt')
FIT_down = np.genfromtxt('td_140.0_down.txt')

if debug1:
    for il in range(len(z_range)):
        bin_width = 0.02
        if pmt == 'up':
            t_offset = z_range[il]*1E-2/(c/n_ref)*1E9
        elif pmt == 'down':
            t_offset = z_range[-1-il]*1E-2/(c/n_ref)*1E9
        bin_offset = round(t_offset/bin_width) 

        xx = np.linspace(0,20,1000)
        y_fit = moffat_king_exponential(xx, par=FIT[il], arg=tau)
        y_fit = np.roll(y_fit, bin_offset)

        fig, ax = plt.subplots(1,1)
        ax.step(T[il], Y[il], where = "mid", color = "k", label = "z = {:.0f} cm".format(z_range[il]))
        ax.plot(xx, y_fit, color = 'blue', ls = '--', 
                label = r'fit, $\chi^2 = {:.2e}$'.format(loss(par = FIT[il], args = [hist_x, fold_y, z_range[il], tau])))
        ax.axvline(T_max[il], color = "purple", label = "maximum")
        ax.axvline(T_50[il], color = "red", label = "median")
        ax.axvspan(T_lhm[il], T_uhm[il], color = "red", label = "FWHM", alpha = 0.5)
        ax.axvspan(T_lq90[il], T_uq90[il], color = "C1", label = "90% quantile", alpha = 0.5)
        ax.axvspan(T_lq95[il], T_uq95[il], color = "C1", label = "95% quantile", alpha = 0.25)
        ax.set_xlabel("propagation time t [ns]")
        ax.set_ylabel("normalized entries")
        ax.set_xlim((0,20))
        ax.set_ylim((0,0.5))
        ax.legend()
        plt.tight_layout()
        plt.savefig("./plots/hist_pmt_"+pmt+"/delay_hist_z={:.0f}cm_".format(z_range[il])+pmt+".png")
        #plt.show()
        plt.close()        

if debug2:
    # sum up all folded histograms
    t = hist_x
    y = np.sum(Y, axis=0)
    y /= np.sum(y)
    # get max, median, quantiles, FWHM of combined distribution
    t_max, t_50, t_lq90, t_uq90, t_lq95, t_uq95, t_lhm, t_uhm = hist_charac(t, y)

    fig, ax = plt.subplots(1,1)
    ax.step(t, y, where = "mid", color = "k", label = "comb. hist.")
    ax.axvline(t_max, color = "purple", label = "maximum")
    ax.axvline(t_50, color = "red", label = "median")
    ax.axvspan(t_lhm, t_uhm, color = "red", label = "FWHM", alpha = 0.5)
    ax.axvspan(t_lq90, t_uq90, color = "C1", label = "90% quantile", alpha = 0.5)
    ax.axvspan(t_lq95, t_uq95, color = "C1", label = "95% quantile", alpha = 0.25)
    ax.set_xlabel("propagation time t [ns]")
    ax.set_ylabel("normalized entries")
    ax.set_xlim((0,20))
    ax.set_ylim((0,0.02))
    ax.legend()
    plt.tight_layout()
    plt.savefig("./plots/delay_hist_comb_"+pmt+".png")
    plt.show()
    plt.close()

if debug3:
    fig, ax = plt.subplots(1,1)
    ax.plot(z_range, T_max, color = "purple", label = "maximum")
    ax.plot(z_range, T_50, color = "red", label = "median")
    ax.fill_between(z_range, T_lhm, T_uhm, color = "red", alpha = 0.5, label = "FWHM")
    ax.fill_between(z_range, T_lq90, T_uq90, color = "C1", alpha = 0.5, label = "90% quantile")
    ax.fill_between(z_range, T_lq95, T_uq95, color = "C1", alpha = 0.25, label = "95% quantile")
    ax.set_xlabel("z position [cm]")
    ax.set_ylabel(r"propagation time [ns]")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig("./plots/time_delay_"+pmt+".png")
    plt.show()

if debug4:
    zz = np.linspace(0,74,200)

    a_ma = moving_average(FIT.T[0], n=5, const_padding = True)
    s_ma = moving_average(FIT.T[1], n=5, const_padding = True)
    g_ma = moving_average(FIT.T[2], n=5, const_padding = True)

    a_int = PchipInterpolator(z_range, FIT.T[0])(zz)
    s_int = PchipInterpolator(z_range, FIT.T[1])(zz)
    g_int = PchipInterpolator(z_range, FIT.T[2])(zz)

    a_maint = PchipInterpolator(z_range, a_ma)(zz)
    s_maint = PchipInterpolator(z_range, s_ma)(zz)
    g_maint = PchipInterpolator(z_range, g_ma)(zz)

    fig, ax = plt.subplots(3,1)

    ax[0].plot(z_range, FIT.T[0], label = 'fit')
    ax[0].plot(z_range, a_ma, label = 'ma', ls = ':')
    ax[0].plot(zz, a_int, label = 'int', ls = '--')
    ax[0].plot(zz, a_maint, label = 'ma-int', ls = '-.')
    ax[0].set_xlabel("z position [cm]")
    ax[0].set_ylabel("Fit parameter A")
    ax[0].set_xlim(0,74)
    ax[0].legend()

    ax[1].plot(z_range, FIT.T[1], label = 'fit')
    ax[1].plot(z_range, s_ma, label = 'ma', ls = ':')
    ax[1].plot(zz, s_int, label = 'int', ls = '--')
    ax[1].plot(zz, s_maint, label = 'ma-int', ls = '-.')
    ax[1].set_xlabel("z position [cm]")
    ax[1].set_ylabel(r"Fit parameter $\sigma$")
    ax[1].set_xlim(0,74)
    ax[1].legend()

    ax[2].plot(z_range, FIT.T[2], label = 'fit')
    ax[2].plot(z_range, g_ma, label = 'ma', ls = ':')
    ax[2].plot(zz, g_int, label = 'int', ls = '--')
    ax[2].plot(zz, g_maint, label = 'ma-int', ls = '-.')
    ax[2].set_xlabel("z position [cm]")
    ax[2].set_ylabel(r"Fit parameter $\gamma$")
    ax[2].set_xlim(0,74)
    ax[2].legend()


    plt.tight_layout()
    #plt.savefig("./plots/fit_params_"+pmt+".png")
    plt.show()

if debug5:
    xx = np.linspace(0,20,1000)
    zz = np.arange(1,75.5,0.5)
    yy = []
    bin_width = 0.02
        
    for iz, z in enumerate(zz):
        t_offset = z_range[il]*1E-2/(c/n_ref)*1E9
        bin_offset = round(t_offset/bin_width) 
    
        # interpolate parameters of Moffat-King fit at position z
        a = moving_average(FIT.T[0], n=5, const_padding = True)
        s = moving_average(FIT.T[1], n=5, const_padding = True)
        g = moving_average(FIT.T[2], n=5, const_padding = True)

        a = PchipInterpolator(z_range, a)(z)
        s = PchipInterpolator(z_range, s)(z)
        g = PchipInterpolator(z_range, g)(z)
        
        # get pdf
        y = moffat_king_exponential(xx, par = (a,s,g), arg = tau)
        y = np.roll(y, bin_offset)
        yy.append(y)

        # plot two neighboring simulated pdfs together with 'new' interpolated pdf
        if iz%1==0 and iz > 1 and iz <= 75:

            fig, ax = plt.subplots(1,1)
            ax.step(xx, yy[iz], where = "mid", color = "k", label = "sim z = {:.1f} cm".format(zz[iz]))
            ax.plot(xx, yy[iz-1], color = 'blue', ls = '--', label = "intp z = {:.1f} cm".format(zz[iz-1]))
            ax.step(xx, yy[iz-2], where = "mid", color = "grey", label = "sim z = {:.1f} cm".format(zz[iz-2]))
            ax.set_xlabel("propagation time t [ns]")
            ax.set_ylabel("normalized entries")
            ax.set_xlim((0,20))
            ax.set_ylim((0,0.5))
            ax.legend()
            plt.tight_layout()
            plt.show()

if debug6:
    zz = np.linspace(0,74,200)
    fig, ax = plt.subplots(3,1)

    ax[0].plot(z_range, FIT_up.T[1], label = r'PMT$_{up}$')
    ax[0].plot(z_range, FIT_down.T[1], label = r'PMT$_{down}$')
    ax[0].plot(z_range, (FIT_up.T[1]+FIT_down.T[1][::-1])/2, label = r'PMT$_{comb}$')
    ax[0].set_xlabel("z position [cm]")
    ax[0].set_ylabel("Fit parameter A")
    ax[0].set_xlim(0,74)
    ax[0].legend()

    ax[1].plot(z_range, FIT_up.T[2], label = r'PMT$_{up}$')
    ax[1].plot(z_range, FIT_down.T[2], label = r'PMT$_{down}$')
    ax[1].plot(z_range, (FIT_up.T[2]+FIT_down.T[2][::-1])/2, label = r'PMT$_{comb}$')
    ax[1].set_xlabel("z position [cm]")
    ax[1].set_ylabel(r"Fit parameter $\sigma$")
    ax[1].set_xlim(0,74)
    ax[1].legend()

    ax[2].plot(z_range, FIT_up.T[3], label = r'PMT$_{up}$')
    ax[2].plot(z_range, FIT_down.T[3], label = r'PMT$_{down}$')
    ax[2].plot(z_range, (FIT_up.T[3]+FIT_down.T[3][::-1])/2, label = r'PMT$_{comb}$')
    ax[2].set_xlabel("z position [cm]")
    ax[2].set_ylabel(r"Fit parameter $\gamma$")
    ax[2].set_xlim(0,74)
    ax[2].legend()


    plt.tight_layout()
    #plt.savefig("./plots/fit_params_"+pmt+".png")
    plt.show()