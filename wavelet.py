import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import ascii
from matplotlib import ticker, cm
from pylab import *
from astroML.time_series import lomb_scargle
import matplotlib.gridspec as gridspec
import peakutils
import matplotlib.ticker as ticker

__author__ = 'Evgeniya Predybaylo'

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${}e^{{{}}}$'.format(a, b)
def lombscargle(x, y, yerr, Pmin=0.5, Pmax=20, res=10000):
    periods = linspace(Pmin, Pmax, res)
    ang_freqs = 2 * pi / periods
    powers = lomb_scargle(x, y, yerr, ang_freqs, generalized=True)
    return periods, powers
def plot_lcper(time1, flux1, periods1, powers1,lctype):
    scatter_kwargs = {"zorder": 100}

    fig4 = plt.figure(num=None, figsize=(7, 4), facecolor='w', dpi=300)

    gs1 = gridspec.GridSpec(2,1)
    gs1.update(left=0.0,right=0.57,hspace=0.0)
    ax3 = plt.subplot(gs1[0, 0])

    gs2 = gridspec.GridSpec(2,1)
    gs2.update(left=0.665, right=0.97, hspace=0.0)
    ax4 = plt.subplot(gs2[0, 0])

    font_size = 'small'


    ax3.set_ylabel('Normalized Flux', fontsize=font_size, style='normal', family='sans-serif')
    ax3.set_xlabel('BJD - 2450000', fontsize=font_size, style='normal', family='sans-serif')
    ax3.set_xlim([min(time1), max(time1)])
    ax3.plot(time1, flux1, color='#000000',lw=0.5)
    ax3.tick_params(axis='both', labelsize=font_size, direction='in', top=True, right=True)

    ax4.set_ylabel('Power', fontsize=font_size, style='normal', family='sans-serif')
    ax4.set_xlabel('Period (d)', fontsize=font_size, style='normal', family='sans-serif')
    mx = max(periods1)
    ax4.set_xlim([0, 20])
    ax4.set_ylim([0, 1.0])
    ax4.plot(periods1, powers1, color='#000000', lw=0.5)
    peaks = peakutils.indexes(powers1, thres=0.15, min_dist=10)
    for el,i in enumerate(peaks):
        ax4.vlines(periods1[peaks], 0, 1.0, color='#000000', lw=0.75, alpha=0.1, linestyles='--')
        # ax2.text(periods1[i] - 0.58, powers1[i] + 0.015, '%5.2f' % periods1[i], horizontalalignment='center', fontsize='x-small', style='normal', family='sans-serif')
        ax4.text(periods1[i], powers1[i] + (el+1)*0.065, '%5.2f' % periods1[i], horizontalalignment='center', fontsize='x-small', style='normal', family='sans-serif')
    ax4.tick_params(axis='both', labelsize=font_size, direction='in', top=True, right=True)

    # -----------------------------

    #tight_layout()
    fig4.savefig('LombScargle_'+lctype+'.pdf', format='pdf', bbox_inches='tight')
    #show()
    close()
def create_sig(savefile_name):
    mint = 0
    maxt = 80.
    Npoints = 3000
    t = linspace(mint, maxt, Npoints, endpoint=False)

    # self.sig1 = np.cos((2 * np.pi/3.0)*self.t)
    sig = zeros_like(t)
    err = zeros_like(t)+1e-10

    amp1 = 1.0
    amp2 = 1.0
    amp3 = 1.0

    use_t = 0
    for rawr in range(len(t)):
        if t[rawr] < 15.:
            use_t = t[rawr]
            T = 1.5
            A = amp1
        if (t[rawr] >= 15.) and (t[rawr] < 45.):
            use_t = t[rawr] - 15.
            T = 6.0
            A = amp2
        if t[rawr] >= 45.:
            use_t = t[rawr] - 45.
            T = 12.
            A = amp3
        #sig[rawr] = A * np.cos((2 * np.pi / T) * use_t)
        sig[rawr] = np.random.normal(A * np.cos((2 * np.pi / T) * use_t), 0.3)

    sig1 = sig

    dat = array([t, sig1, err])
    dat = dat.T
    fileid = open(savefile_name, 'w+')
    savetxt(fileid, dat, fmt=['%0.6f\t', '%0.6f\t', '%0.6f'])
    fileid.close()

    return sig1
def create_sig_overlap(savefile_name):
    mint = 0
    maxt = 80.
    Npoints = 3000
    t = linspace(mint, maxt, Npoints, endpoint=False)

    # self.sig1 = np.cos((2 * np.pi/3.0)*self.t)
    sig1 = zeros_like(t)
    sig2 = zeros_like(t)
    sig3 = zeros_like(t)
    err = zeros_like(t)+1e-10

    amp1 = 1.0
    amp2 = 1.0
    amp3 = 1.0

    use_t = 0
    for rawr in range(len(t)):
        if t[rawr] < 30.:
            use_t = t[rawr]
            T = 1.5
            A = amp1
            sig1[rawr] = A * np.cos((2 * np.pi / T) * use_t)
            # sig1[rawr] = np.random.normal(A * np.cos((2 * np.pi / T) * use_t), 0.3)
        if (t[rawr] >= 10.) and (t[rawr] < 60.):
            use_t = t[rawr] - 10.
            T = 6.0
            A = amp2
            sig2[rawr] = A * np.cos((2 * np.pi / T) * use_t)
            # sig2[rawr] = np.random.normal(A * np.cos((2 * np.pi / T) * use_t), 0.3)
        if t[rawr] >= 40.:
            use_t = t[rawr] - 40.
            T = 12.
            A = amp3
            sig3[rawr] = A * np.cos((2 * np.pi / T) * use_t)
            # sig3[rawr] = np.random.normal(A * np.cos((2 * np.pi / T) * use_t), 0.3)

        sig = sig1 + sig2 + sig3

    sig1 = sig

    dat = array([t, sig1, err])
    dat = dat.T
    fileid = open(savefile_name, 'w+')
    savetxt(fileid, dat, fmt=['%0.6f\t', '%0.6f\t', '%0.6f'])
    fileid.close()

    return sig1
def read_sig_data(readfile_name):
    filedata = ascii.read(readfile_name, data_start=0)
    time = np.array(filedata['col1'])
    data = np.array(filedata['col2'])
    err = np.array(filedata['col3']) + 1e-10

    # ------------ Calculate Variance ----------------
    data = data - np.mean(data)
    variance = np.std(data, ddof=1) ** 2
    print("variance = ", variance)

    return time, data, err, variance
def read_model_data(readfile_model_name,planet=True):
    filedata = ascii.read(readfile_model_name, data_start=0)
    time = np.array(filedata['col1'])
    if planet == True:
        data = np.array(filedata['col2'])
    if planet == False:
        data = np.array(filedata['col3'])
    err = np.zeros_like(data) + 1e-10

    # ------------ Calculate Variance ----------------
    data_for_var = data - np.mean(data)
    variance = np.std(data_for_var, ddof=1) ** 2
    print("variance = ", variance)

    return time, data, err, variance
def determine_levels(power):

    levels = np.array([0., 2., 4., 8., 16., 32., 64., 128., 256.])
    maxpow = np.max(power)
    if (maxpow <= 100) and (maxpow > 10):
        levels = levels * 1e-1
    if (maxpow <= 10) and (maxpow > 1):
        levels = levels * 1e-2
    if (maxpow <= 1) and (maxpow > 0.1):
        levels = levels * 1e-3
    if (maxpow <= 0.1) and (maxpow > 0.01):
        levels = levels * 1e-4
    if (maxpow <= 0.01) and (maxpow > 0.001):
        levels = levels * 1e-5
    if (maxpow <= 0.001) and (maxpow > 0.0001):
        levels = levels * 1e-6

    # fix levels according to max power
    while maxpow < levels[-1]:
        levels = levels[0:-1]

    return levels
def plot_wavelet_power_spectrum(time,data,period,power,global_ws,lctype,wavelet_type):

    #----------------------------------------------------------------------#
    #----------------------------Plotting----------------------------------#
    #----------------------------------------------------------------------#

    font_size = 'small'

    levels = determine_levels(power)

    #--- Plot time series
    fig = plt.figure(figsize=(6.75, 5.5), dpi=300)
    gs = GridSpec(2, 4, hspace=0.30, wspace=0.10)
    plt.subplots_adjust(left=0.1, bottom=0.08, right=0.95, top=0.95, wspace=0, hspace=0)
    plt.subplot(gs[0, 0:4])
    plt.plot(time, data, 'k', lw=1.0)
    plt.xlim(xlim[:])
    plt.xlabel('Time', fontsize=font_size)
    plt.ylabel('Data', fontsize=font_size)
    plt.title('(a) Data', fontsize='x-small')
    plt.tick_params(axis='both',labelsize=font_size, direction='in', top=True, right=True)

    # plt.text(time[-1] + 35, 0.5,'Wavelet Analysis\nC. Torrence & G.P. Compo\n' +
    #     'http://paos.colorado.edu/\nresearch/wavelets/',
    #     horizontalalignment='center', verticalalignment='center')



    #--- Contour plot wavelet power spectrum
    # plt3 = plt.subplot(3, 1, 2)
    plt3 = plt.subplot(gs[1, 0:3])
    #levels = np.concatenate((np.arange(0,np.max(power),2),[999.]))
    # levels = arange(0,np.max(power),1)
    #levels = np.arange(0,np.max(power),50)
    # CS = plt.contourf(time, period, power, len(levels))  #*** or use 'contour'
    # im = plt.contourf(CS, levels=levels, colors=['white','bisque','orange','orangered','darkred'])

    time_plot,period_plot = np.meshgrid(time, period)
    boundaries = levels
    cmap_reds = plt.cm.get_cmap('Reds',len(boundaries)) # create list of colors from colormap with length of levels array - 1
    colors = list(cmap_reds(np.arange(len(boundaries))))
    colors[0] = "white" #replace first color with white
    cmap = matplotlib.colors.ListedColormap(colors[:-1], "")
    cmap.set_over(colors[-1]) # set over-color to last color of list
    norm = mpl.colors.BoundaryNorm(boundaries, ncolors=len(boundaries)-1, clip=False)
    CS = plt.pcolormesh(time_plot, period_plot, power, norm=norm, cmap=cmap, vmax=np.max(levels))
    fig.colorbar(CS, extend='max',format='%.0e') # format=ticker.FuncFormatter(fmt))
    plt.xlabel('Time', fontsize=font_size)
    plt.ylabel('Period', fontsize=font_size)
    plt.title('(b) Wavelet Power Spectrum', fontsize='x-small')
    plt.xlim(xlim[:])
    # 95# significance contour, levels at -99 (fake) and 1 (95# signif)
    plt.contour(time, period, sig95, [-99, 1], colors='k')
    # cone-of-influence, anything "below" is dubious
    plt.plot(time, coi, 'k')
    # --- format y-scale ---
    plt3.set_yscale('log', basey=2, subsy=None)
    plt.ylim([np.min(period), np.max(period)])
    ax = plt.gca().yaxis
    ax.set_major_formatter(ticker.ScalarFormatter())
    #plt3.ticklabel_format(axis='y', style='plain')
    plt3.invert_yaxis()
    plt3.tick_params(axis='both',labelsize=font_size, direction='in', top=True, right=True)
    # set up the size and location of the colorbar
    # position=fig.add_axes([0.5,0.36,0.2,0.01])
    # plt.colorbar(im, cax=position, orientation='horizontal') #, fraction=0.05, pad=0.5)

    # plt.subplots_adjust(right=0.7, top=0.9)

    #--- Plot global wavelet spectrum
    plt4 = plt.subplot(gs[1, -1])
    plt.plot(global_ws, period, lw=1.0, color='#000000')
    signif_color = colors[-2]
    plt.plot(global_signif, period, '--', lw=1.0, color=signif_color)
    peaks = peakutils.indexes(global_ws, thres=0.05, min_dist=10)
    for el, i in enumerate(peaks):
        plt4.hlines(period[peaks], 0, 1.25 * np.max(global_ws), color='#000000', lw=0.75, alpha=0.1, linestyles='--')
        plt4.text( global_ws[i] + ((1./256.)*global_ws[i] * 4), period[i], '%5.2f' % period[i], horizontalalignment='center',fontsize='x-small', style='normal', family='sans-serif', rotation=270)
    plt.xlabel('Power', fontsize=font_size)
    plt.title('(c) Global Wavelet Spectrum', fontsize='x-small')
    plt.xlim([0, 1.35 * np.max(global_ws)])
    # format y-scale
    plt4.set_yscale('log', basey=2, subsy=None)
    plt.ylim([np.min(period), np.max(period)])
    ax = plt.gca().yaxis
    ax.set_major_formatter(ticker.ScalarFormatter())
    #plt4.ticklabel_format(axis='y', style='plain')
    plt4.invert_yaxis()
    plt4.tick_params(axis='both',labelsize=font_size, direction='in', top=True, right=True)
    plt4.set_yticklabels([])

    # --- Plot 2--8 yr scale-average time series
    # plt.subplot(gs[2, 0:3])
    # plt.plot(time, scale_avg, 'k')
    # plt.xlim(xlim[:])
    # plt.xlabel('Time', fontsize=font_size)
    # plt.ylabel('Avg variance', fontsize=font_size)
    # plt.title('(d) 2-8 yr Scale-average Time Series', fontsize=font_size)
    # plt.plot(xlim, scaleavg_signif + [0, 0], '--')
    fig.savefig('wavelet_'+lctype+'_'+wavelet_type+'.png', format='png') #, bbox_inches='tight')
    #plt.show()
    close()
def plot_loglog_gradient(period, global_ws, lctype, wavelet_type, star_rot, planet_orb):

    fig4 = plt.figure(num=None, figsize=(7, 6), facecolor='w', dpi=300)

    font_size = 'small'

    plt.title(wavelet_type, fontsize=font_size, style='normal', family='sans-serif')

    gs1 = gridspec.GridSpec(2,1)
    gs1.update(left=0.0,right=0.44,hspace=0.0)
    ax3 = plt.subplot(gs1[0, 0])

    gs2 = gridspec.GridSpec(2,1)
    gs2.update(left=0.53, right=0.97, hspace=0.0)
    ax4 = plt.subplot(gs2[0, 0])

    # convert period to microHz
    freq_days = 1.0/period
    Hz = freq_days*(1./86400.)
    microHz = Hz * 1e-6

    delta_freq = np.diff(microHz)
    mid_freq = microHz[0:-1] + 0.5*delta_freq

    norm_global_ws = global_ws / np.max(global_ws)

    global_ws_kplus1 = norm_global_ws[1:len(norm_global_ws)]
    global_ws_k = norm_global_ws[0:-1]

    dlnP = np.diff(np.log(norm_global_ws))
    dlnfreq = np.diff(np.log(microHz))

    grad = 1. + (dlnP/dlnfreq) * (delta_freq)/mid_freq
    grad2 = global_ws_kplus1/global_ws_k

    star_rot_freq = (1.0 / star_rot) * (1./86400.) * 1e-6
    planet_orb_freq = (1.0 / planet_orb) * (1. / 86400.) * 1e-6


    ax3.set_ylabel('Normalized Power', fontsize=font_size, style='normal', family='sans-serif')
    ax3.set_xlabel(r'Frequency ($\mu$Hz)', fontsize=font_size, style='normal', family='sans-serif')
    ax3.set_xlim([np.min(microHz), np.max(microHz)])
    ymin = 1e-6
    ymax = 1.25*np.max(norm_global_ws)
    ax3.set_ylim([ymin, ymax])
    ax3.plot(microHz, norm_global_ws, color='blue',lw=1.0, label=wavelet_type)

    ax3.vlines([star_rot_freq], ymin, ymax, color='#000000', lw=0.75, linestyles='--')
    ax3.vlines([planet_orb_freq], ymin, ymax, color='#000000', lw=0.75, linestyles='--')
    ax3.text(star_rot_freq, 2e-6, '%5.1f' % star_rot, horizontalalignment='left', fontsize='x-small',style='normal', family='sans-serif')
    ax3.text(planet_orb_freq, 2e-6, '%5.1f' % planet_orb, horizontalalignment='left', fontsize='x-small',style='normal', family='sans-serif')

    ax3.tick_params(axis='both', labelsize=font_size, direction='in', top=True, right=True)
    ax3.set_yscale('log', basey=10, subsy=None)
    ax3.set_xscale('log', basex=10, subsx=None)
    ax3.legend(loc='upper right',fontsize=font_size)


    ax4.set_ylabel('Gradient', fontsize=font_size, style='normal', family='sans-serif')
    ax4.set_xlabel(r'Frequency ($\mu$Hz)', fontsize=font_size, style='normal', family='sans-serif')
    ax4.set_xlim([np.min(microHz), np.max(microHz)])
    yminlim = np.min(grad)-0.15
    ymaxlim = np.max(grad)+0.15
    ax4.set_ylim([yminlim, ymaxlim])
    # ax4.plot(microHz, grad, color='red', lw=0.5)
    ax4.plot(mid_freq, grad, color='#ffa31a', lw=1.0, label='Method 1')
    ax4.plot(mid_freq, grad2, color='#2eb82e', lw=1.0, label='Method 2')
    peaks = peakutils.indexes(grad, thres=0.15, min_dist=10)
    for el,i in enumerate(peaks):
        #import pdb; pdb.set_trace()
        ax4.vlines(mid_freq[peaks], yminlim, ymaxlim, color='#000000', lw=0.75, alpha=0.1, linestyles='--')
        # ax2.text(periods1[i] - 0.58, powers1[i] + 0.015, '%5.2f' % periods1[i], horizontalalignment='center', fontsize='x-small', style='normal', family='sans-serif')
        peak_per = np.float(mid_freq[i])*1e6*86400.
        ax4.text(mid_freq[i], np.min(grad)-0.125, '%5.2f' % peak_per, horizontalalignment='left', fontsize='x-small', style='normal', family='sans-serif')
    ax4.tick_params(axis='both', labelsize=font_size, direction='in', top=True, right=True)
    ax4.set_xscale('log', basex=10, subsx=None)
    ax4.legend(loc='upper right',fontsize=font_size)
    # -----------------------------

    #tight_layout()
    fig4.savefig('log-log-gradient_'+lctype+'_'+wavelet_type+'.pdf', format='pdf', bbox_inches='tight')
    close()


# WAVETEST Example Python script for WAVELET, using NINO3 data dataset
#
# See "http://paos.colorado.edu/research/wavelets/"
# The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
#
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
#   changed all "log" to "log2", changed logarithmic axis on GWS to
#   a normal axis.
# ------------------------------------------------------------------------------------------------------------------

# READ THE DATA
# data = np.loadtxt('sst_nino3.dat')  # input data time series


readfile_sig_name = 'test_sig_overlap.txt'
readfile_model_name = 'model_lightcurves/lc_4.txt'

lctype = 'model' # 'sig'
wavelet_type = 'PAUL'


overlap = True
scatter = True
if lctype == 'sig':
    # Need to implement easy scatter and overlap and amplitude test options
    create_sig('test_sig_scatter.txt')
    create_sig_overlap('test_sig_overlap.txt')
    time, data, err, variance = read_sig_data(readfile_sig_name)
    filesave_figname = 'wavelet_test_sig_overlap.png'
if lctype == 'model':
    filesave_figname = 'wavelet_test_model.png'
    filesave_log_figname = 'log-log-gradient_model.pdf'
    time, data, err, variance = read_model_data(readfile_model_name,planet=False)
    stellar_rotation = 7 # days
    planet_orbit = 12 # days


#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
if 0:
    variance = 1.0
    data = data / np.std(data, ddof=1)
n = len(data)
# dt = 0.25
# time = np.arange(len(data)) * dt + 1871.0  # construct time array
# xlim = ([1870, 2000])  # plotting range
dt = np.mean(np.diff(time))
xlim = ([np.min(time),np.max(time)])
pad = 1  # pad the time series with zeroes (recommended)
#dj = 0.25  # this will do 4 sub-octaves per octave
s0 = 2 * dt  # this says start at a scale of twice the time difference between data points
#j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.1  # lag-1 autocorrelation for red noise background 0.72
print("lag1 = ", lag1)
mother = wavelet_type


#------------ LOMB-SCARGLE APPROACH ----------------
pers, pows = lombscargle(time, data, err, Pmin=s0, Pmax=np.max(time)-np.min(time))
plot_lcper(time, data, pers, pows,lctype=lctype)
peakpows_LS = peakutils.indexes(pows, thres=0.2, min_dist=10)
peakpers_LS = pers[peakpows_LS]


#------------ WAVELET TRANSFORM ----------------
# wave, period, scale, coi = wavelet(data, dt, pad, dj, s0, j1, mother='MORLET')
wave, period, scale, coi = wavelet(data, dt, pad, dj=0.05 ,mother=mother,param=6)
power = (np.abs(wave)) ** 2  # compute wavelet power spectrum
global_ws = (np.sum(power, axis=1) / n)  # time-average over all times

# Significance levels:
signif = wave_signif(([variance]), dt=dt, sigtest=0, scale=scale,
    lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
sig95 = power / sig95  # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1,
    lag1=lag1, dof=dof, mother=mother)

# Scale-average between El Nino periods of 2--8 years
# avg = np.logical_and(scale >= 2, scale < 8)
# Cdelta = 0.776  # this is for the MORLET wavelet
# scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
# scale_avg = power / scale_avg  # [Eqn(24)]
# scale_avg = dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
# scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2,
#     lag1=lag1, dof=([2, 7.9]), mother=mother)

plot_wavelet_power_spectrum(time,data,period,power,global_ws,lctype,wavelet_type)

plot_loglog_gradient(period, global_ws, lctype, wavelet_type,star_rot=stellar_rotation,planet_orb=planet_orbit)





