import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.timeseries import LombScargle
import scipy.signal as signal
import numpy as np
import pandas as pd

data = pd.read_csv('2017-18.dat', sep = '\s+', header = None)

# search and remember NaN in data into logical massive mask
maska = (np.isnan(data[1].astype(float)))
for i in range(len(maska)): 
    maska[i] = (not maska[i])

# skip NaN in data
time = np.asarray(data[maska][0])
foEs = np.asarray(data[maska][1])
ntime = len(time)

# Процедура ко всему ряду
foEsmean = np.nanmean(foEs)

## ------------ Считаем по всему ряду -----
# min and max frequensy in 1/day
min_freq = 1.0 / time[ntime - 1]
max_freq = 1.0
freq, Y = LombScargle(time, foEs - foEsmean).autopower(minimum_frequency=min_freq, maximum_frequency=max_freq)
fig, ax = plt.subplots(2, 1)
ax[0].plot(time, foEs)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')

ax[1].plot(freq, Y, 'r') # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('Y(freq)')

plt.show()

## -----------------------------------------

plt.xlim([2.,30.])
plt.plot(1.0 / freq, Y, 'r')
plt.xlabel('Period(day)')
plt.ylabel('Y')

plt.show()

# готовим оконное преобразование
# window size in day
w_size = 18.0

dx = 0.042 # step into win
# number of points in windows

w_num = int(w_size / dx)
# number of time's points
w_frec = len(time) - 1

# возможны пропуски в данных, ищем их (вроде после удаленгия NaN не надо, но осталось)
bord=[]
i_st_end = np.zeros(2).astype(int)
jstart = 0
jend = 0

# search windows
for i in range(0, ntime - w_num - 1): # цикл по всем точкам xt
    ix_start = i
    ix_end = i + w_num
    x_s = time[ix_start]
    x_e = time[ix_end]
    while (x_s > time[jstart]):
        jstart = jstart+1
        jend = jstart
    while (x_e > time[jend]):
        jend = jend + 1
        if (jend > w_frec):
            jend = w_frec
            break

    bord.append([jstart,jend])

w_spectr = len(bord)

minW_freq = 1.0 / w_size
maxW_freq = 1.0 / 0.05
freq0, Y = LombScargle(time[0:w_num], foEs[0:w_num]).autopower(minimum_frequency=minW_freq, maximum_frequency=maxW_freq)

i_df = len(freq0)
m_LS=np.zeros(i_df * w_spectr).reshape(i_df, w_spectr)

# massive of spectrogramm
# window LS transform
i = 0
for gran in bord:
    i_st=gran[0]
    i_en=gran[1]
    foEsmean=np.nanmean(foEs[i_st:i_en])
    
    m_LS[0:i_df,i] = LombScargle(time[i_st:i_en], foEs[i_st:i_en] - foEsmean).power(freq0)
    i += 1

plt.contourf(time[1:-w_num-1],1./freq0,m_LS[:,1:], levels=20, cmap=cm.rainbow)

plt.colorbar()
plt.show()