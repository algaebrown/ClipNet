import numpy as np
import matplotlib.pyplot as plt
def filter_key_transcript(key_transcript, read_density):
    key_transcript = key_transcript.merge(c = '3,6,7', o = 'distinct', s=True)
    good_signal = []
    for i in range(len(key_transcript)):
        chro = key_transcript[i].chrom
        start = int(key_transcript[i].start)
        end = int(key_transcript[i].stop)
        strand = key_transcript[i].strand
        if np.nansum(read_density.values(chro, start, end, strand)) < 0.01:
            pass
        else:
            good_signal.append(key_transcript[i])
    return good_signal
        
    
        
            
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

   
    x = np.nan_to_num(x,0)
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y
def plot_waveform(key_transcript, read_density, is_smooth = False, window_len = None, sharex = True):
    f, ax = plt.subplots(4,10, figsize = (20,16), sharex = sharex)
    ax = ax.flatten()
    for i in range(40):
        chro = key_transcript[i].chrom
        start = int(key_transcript[i].start)
        end = int(key_transcript[i].stop)
        strand = key_transcript[i].strand
    
        shape = read_density.values(chro, start, end, strand)
        if strand == '-':
            shape = [-s for s in shape]
        if is_smooth:
            ax[i].plot(smooth(shape, window_len = window_len))
        else:
            ax[i].plot(shape)
        ax[i].set_title('{0}:{1}-{2}, {3}'.format(chro, start, end, strand))

from scipy.fft import rfft
def plot_spectral(key_transcript, read_density):
    f, ax = plt.subplots(10,1, figsize = (2,10), sharex = True)
    
    for i in range(10):
        chro = key_transcript[i].chrom
        start = int(key_transcript[i].start)
        end = int(key_transcript[i].stop)
        strand = key_transcript[i].strand
        
        shape = read_density.values(chro, start, end, strand)
        if strand == '-':
            shape = [-s for s in shape]
        shapef = rfft(np.nan_to_num(shape,0)/np.nansum(shape))
        ax[i].plot(np.abs(shapef))
        ax[i].set_yscale('log')
        ax[i].set_title('{0}:{1}-{2}, {3}'.format(chro, start, end, strand))
import math
def spectrogram(shape, step = 'constant', step_size = 1000):
    '''
    step = constant, equal step size; relative, cut into 10 segments
    '''
    shape = np.nan_to_num(shape,0)
    # normalize
    shape = shape/np.sum(shape)
    
    if step == 'constant':
        walk = int(step_size/2)
        
        
    if step == 'relative':
        step_size = math.floor(len(shape)/8)
        walk = math.floor(step_size/2)
    i,j = 0, step_size
    window = np.hamming(step_size)
    all_timepoint = []
    while j < len(shape):
        
        segment = shape[i:j]
        y = window*segment
        segf = rfft(y)
        all_timepoint.append(np.abs(segf))
        i+= walk
        j += walk
        
    return np.stack(all_timepoint)
import seaborn as sns
def plot_spectrogram(key_transcript, read_density, step = 'relative'):
    f, ax = plt.subplots(2,5, figsize = (20,10), sharex = True)
    ax = ax.flatten()
    
    for i in range(10):
        chro = key_transcript[i].chrom
        start = int(key_transcript[i].start)
        end = int(key_transcript[i].stop)
        strand = key_transcript[i].strand
    
        shape = read_density.values(chro, start, end, strand)
        print(chro, start, end, strand)
        if strand == '-':
            shape = [-s for s in shape]
        spect = spectrogram(shape, step = step)
        sns.heatmap(spect[:, :50], yticklabels = np.arange(start = start, stop = end, step = math.floor((end-start)/16)), vmax = 0.2, ax = ax[i], cbar = False)
        ax[i].set_title('{0}:{1}-{2}, {3}'.format(chro, start, end, strand))