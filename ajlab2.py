import numpy as np
import os
from scipy import signal

def calc_coh_from_dat (dat_filename,num_channels,sample_rate,t0=0,chuck_size_in_seconds=10,freqs=[300,500],filt=False):
    #map = np.concatenate((np.array([29,19,18,28,30,20,17,21,31,22,16,23,27,26,25,24,7,6,5,4,8,10,9,3,11,2,12,1,13,0,14,15]), 
    #                       np.array([29,19,18,28,30,20,17,21,31,22,16,23,27,26,25,24,7,6,5,4,8,10,9,3,11,2,12,1,13,0,14,15]) + 32))
    map=np.concatenate((np.array([20,19,25,18,24,17,22,32,21,31,26,8,23,13,30,12,29,11,28,9,27,3,4,15,7,14,10,2,6,1,5,16]),np.array([20,19,25,18,24,17,22,32,21,31,26,8,23,13,30,12,29,11,28,9,27,3,4,15,7,14,10,2,6,1,5,16])+32))
    #map=np.concatenate((np.array([13,14,8,15,9,16,11,1,12,2,7,25,10,20,3,21,4,22,5,24,6,30,29,18,26,19,23,31,27,32,28,17]),np.array([13,14,8,15,9,16,11,1,12,2,7,25,10,20,3,21,4,22,5,24,6,30,29,18,26,19,23,31,27,32,28,17])+32))
    map-=1

    #map=np.concatenate((np.array([29,27,25,23,19,18,17,16,20,22,21,24,26,28,30,31,10,6,2,1,5,7,9,11,15,14,13,12,3,0,4,8]),np.array([29,27,25,23,19,18,17,16,20,22,21,24,26,28,30,31,10,6,2,1,5,7,9,11,15,14,13,12,3,0,4,8])+32))
    #map=np.concatenate((np.array([30,28,26,24,20,19,18,17,21,23,22,25,27,29,31,32,11,7,3,2,6,8,10,12]),np.array([30,28,26,24,20,19,18,17,21,23,22,25,27,29,31,32,11,7,3,2,6,8,10,12])+32))
    #map=np.concatenate((np.array([32,31,29,27,25,22,23,21,17,18,19,20,24,26,28,30,9,5,1,4,13,14,15,16,12,10,8,6,2,3,7,11]),np.array([32,31,29,27,25,22,23,21,17,18,19,20,24,26,28,30,9,5,1,4,13,14,15,16,12,10,8,6,2,3,7,11])+32))
    #map=np.concatenate((np.array([9,5,1,4,13,14,15,16,12,10,8,6,2,3,7,11,32,31,29,27,25,22,23,21,17,18,19,20,24,26,28,30]),np.array([9,5,1,4,13,14,15,16,12,10,8,6,2,3,7,11,32,31,29,27,25,22,23,21,17,18,19,20,24,26,28,30])+32))
    if (filt):
        nyq = 0.5 * sample_rate
        low = freqs[0] / nyq
        high = freqs[1]/ nyq    
        filtpar_b, filtpar_a = signal.butter(3,  [low, high], btype='band')

    statinfo = os.stat(dat_filename)
    filesize=statinfo.st_size
    pop_sample_size_in_bytes = num_channels*2
    num_samples = int(filesize / pop_sample_size_in_bytes)
    mem_map_data = np.memmap(dat_filename, dtype='int16', mode='r', shape=(num_samples, num_channels))
    mem_map_data.shape
    rec_time = mem_map_data.shape[0]/sample_rate
    #print (rec_time)
    i0=t0*sample_rate
    i_width=chuck_size_in_seconds*sample_rate
    i1=i0+i_width
    data = 10e-6 * 0.195 * mem_map_data[int(i0):int(i1), :]
    data = data[:, map]


    count = 0
    subsampling=20
    count = 0;
    C=np.zeros(shape=(num_channels,num_channels), dtype='float')
    print ("Subsampling %d " % subsampling)

    for ch1 in range(0,num_channels):
        for ch2 in range(ch1+1,num_channels):
            x=data[:,ch1]
            if (filt):
                y = signal.lfilter(filtpar_b, filtpar_a, x)
                x1 = y[::subsampling]
            else:
                x1 = x[::subsampling]
            x=data[:,ch2]
            if (filt):
                y = signal.lfilter(filtpar_b, filtpar_a, x)
                x2 = y[::subsampling]
            else:
                x2 = x[::subsampling]
            f, Cxy = signal.coherence(x1, x2, int(sample_rate/subsampling), nperseg=1024)
            C[ch1,ch2]=np.mean(Cxy[(f >= freqs[0]) & (f <= freqs[1])])

    del mem_map_data
    return (C)
