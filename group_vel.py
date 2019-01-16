# -*- coding: utf-8 -*-
# imports
import glob
from colorama import init, Style, Fore
import math
#init(autoreset=True)
from matplotlib import pyplot as plt
import numpy as np
np.seterr(divide='ignore')


# print(Style.BRIGHT + Fore.RED + 'Oops, invalid Syntax. Try again with a comma maybe?!')




"""

# DO:
# ---
# Calcul the dispersion curve using a gaussian filter
# and correct of the freq using centroid freq
# plot it
# Param:
# ------
"""


def multi_filter(frequencies, cross_spectrum, filter_freqs, bandwidth, nzeropad=0):
    """
    Function that filter the signal with a set of gaussian windows.
    It applies the centroid correction from Campillo et al., 1996.
    Give you the main input useful for the FTAN and then picking
    of the Group velocoties.

    :type frequencies: :class:`~numpy.ndarray`
    :param frequencies: 1-D array containing frequency samples of the CC spectrum.
    :type corr_spectrum: :class: `~numpy.ndarray`
    :param corr_spectrum: 1-D or 2-D array containing real or complex CC spectrum.
    :type filter_freqs: :class:`~numpy.ndarray`
    :param filter_freqs: 1-D array containing the target filter for FTAN.
    :type bandwidth: :class: float
    :param bandwidth: parameter (sigma) defining the width of the gaussian filter.
    :type nzeropad
    :param nzeropad

    :rtype: :class:`~numpy.ndarray`
    :return: **centroid_freqs, signal_tf** - centroid_freqs is the corrected frequency
        set and signal_tf is the filtered signal at each centroid_freqs.
    """

    class Gauss:
        def __init__(self, f0, a=1.0):
            self._omega0 = 2. * math.pi * f0
            self._a = a

        def evaluate(self, freqs):
            omega = 2. * math.pi * freqs
            return np.exp(-((omega - self._omega0)
                   / (self._a * self._omega0)) ** 2)[0]

    ninit = len(cross_spectrum)
    ntot = 2 * (ninit - 1) + nzeropad
    nfilt = len(filter_freqs)
    if np.shape(frequencies) != (1, ninit):
        frequencies = np.expand_dims(frequencies, axis=0)

    signal_tf = np.zeros((nfilt, ntot))
    centroid_freqs = np.zeros(nfilt)
    for ifilt, f0 in enumerate(filter_freqs):
        taper = Gauss(f0, a=bandwidth)
        weights = taper.evaluate(frequencies)
        analytic_spec = np.zeros(ntot, dtype=np.complex)
        analytic_spec[:ninit] = cross_spectrum * weights

        enorm = np.abs(analytic_spec[:ninit]) ** 2
        enorm /= np.sum(enorm)

        if ninit % 2 == 0:
            analytic_spec[1:ninit - 1] *= 2.
        else:
            analytic_spec[1:ninit] *= 2.

        analytic = np.fft.ifft(analytic_spec)
        signal_tf[ifilt, :] = np.abs(analytic)

        enorm = np.abs(analytic_spec[:ninit]) ** 2
        enorm /= np.sum(enorm)
        centroid_freqs[ifilt] = np.sum(frequencies * enorm)

    return centroid_freqs, signal_tf


def ftan_vg(frequencies, cross_spectrum, interstation_distance, umin, umax,
            freqmin, freqmax, sigma, nfilt, nzeropad=0, normalize_amp=True,
            plot=True):

    # original values
    nf = len(frequencies) # original freqs
    df = frequencies[1] - frequencies[0] #original df

    # new values with nzeropad
    nt = 2 * (nf - 1) + nzeropad # new n sample
    ndf = 1 / ((2 * (nf - 1) + nzeropad) * df) # new df
    ntvec = np.arange(0,  ndf * (2 * (nf - 1) + nzeropad), ndf) # new time vec

    # new velgroups axe
    velgroups = interstation_distance / ntvec
    nt2 = np.min([math.floor(nt / 2) + 1, np.where(velgroups > umin)[0][-1]])# to cut signal for plot
    nt1 = np.where(velgroups < umax)[0][0] # to cut signal for plot
    velgroups = velgroups[nt1:nt2]

    filter_freqs = np.linspace(freqmin, freqmax, nfilt)
    centroid_freqs, signal_tf = multi_filter(frequencies, cross_spectrum, filter_freqs, sigma, nzeropad)

    # cutting the signal
    signal_tf = signal_tf[:, nt1:nt2]

    # calcul of the dispersion curve = Vg for the max of amplitude at each c_freq
    vg_picks = []
    for tf in signal_tf:
        if normalize_amp:
            tf /= np.max(tf)

        vg_picks.append(velgroups[np.argmax(tf)])

    print(np.shape(vg_picks),np.shape(centroid_freqs))




    if plot:
        Z = np.transpose(signal_tf)
        X, Y = np.meshgrid(centroid_freqs, velgroups)

        fig = plt.figure(figsize=(6,4))
        plt.title('FTAN Spectrum')
        plt.ylabel('$U$ [$km/s$]')
        plt.xlabel('Freq [$Hz$]')
        plt.ylim(umin, umax)
        plt.xlim(centroid_freqs.min(), centroid_freqs.max())
        fac0 = 1e-3
        #step = max(round(nt2/600),1)
        Z = np.log((fac0 * np.abs(Z)) + 1)
        plt.contourf(X, Y, Z,256, cmap=plt.jet())
        plt.scatter(centroid_freqs, vg_picks)
        plt.show()



def greens(tr, peg=True):
    # Cutt the trace
    tr.append(tr.ydata[0])
    left = tr.chop(tr.tmin, tr.tmin + ((len(tr.ydata) / 2) * tr.deltat), inplace=False, include_last=True)
    left.set_ydata(left.get_ydata()[::-1])  # reverse data
    left.shift(-tr.tmin)  # let it begin at t=0
    right = tr.chop(tr.tmin + ((len(tr.ydata) / 2) * tr.deltat), tr.tmax, inplace=False, include_last=True)
    if peg:
        plug = tr.copy()
        # print plug, plug.ydata
        # print len(right.ydata), len(left.ydata)
        plug.set_ydata(right.ydata + left.ydata[1:])
        tr.set_ydata(plug.ydata)
        # plug.snuffle()
    else:
        if abs(max(left.ydata)) > abs(max(right.ydata)):
            # station = 'LGG%s' %tr.station[-2:]
            # event = 'LGG%s' %tr.station[:2]
            tr = left
        else:
            # station = 'LGG%s' %tr.station[:2]
            # event = 'LGG%s' %tr.station[-2:]
            tr = right
    # for sta in stations:
    #    if station == sta.station:
    #        station = sta
    #    if event == sta.station:
    #        event = sta
    #    tr.plot()
    # dist = orthodrome.distance_accurate50m(event, station)
    dist = 160
    dist /= 1000.
    #    print  Style.BRIGHT + Fore.GREEN + '..\n...\n>>> Distance between station %s and %s = %g km' % (station.station, event.station, round(dist,2))
    tspan = 0
    #   tr.set_station('%s-%s'%(station.station[3:],event.station[3:]))
    # filename=station.station+event.station

    return tr, dist, tspan  # , filename




def main(file, ax):
    umin, umax, Tmin, Tmax, alpha = .001, 1., 1, 25, 0.005  # ask_default()

    # stations = model.load_stations('stations.txt')

    # names =  os.path.basename(file).split('.')
    # names = names[0].split('_')
    station = '001'  # names[1]
    event = '002'  # names[3]
    # event = event.replace('.','')
    # station = station.replace('.','')
    filename = station + '_' + event
    tr = io.load(file, format='from_extension')[0]
    tr, dist, tspan = greens(tr, peg=True)
    print(tr, dist, tspan, umin, umax, Tmin, Tmax, alpha, station, event, filename)
    dispersion(ax, tr, dist, tspan, umin, umax, Tmin, Tmax, alpha, station, event, filename)


if __name__ == '__main__':
    #fig = plt.figure()
    #ax = plt.subplot(111)
   #dd file = '/Users/tzompantli/Desktop/test.mseed'  # /media/VOLC_DATA/Colima/STACKS/01/REF/tmp/AZ_BZN_CI_MPP.mseed'
    #main(file, ax)
    from scipy.io import loadmat as lm
    x = lm('cc.mat')

    frequencies = x['freq']
    frequencies = np.squeeze(x['freq'])
    cross_spectrum = np.squeeze(x['smoothed'])
    interstation_distance=0.16
    alpha = 0.25
    umin =0.1
    umax =1
    freqmin = 1
    freqmax = 25
    nfilt= 100
    nzeropad = 20000

    ftan_vg(frequencies, cross_spectrum, interstation_distance, umin, umax,
            freqmin, freqmax, alpha, nfilt, nzeropad)

    plt.show()

