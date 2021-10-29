from scipy.integrate import simps
from stella_data import time, ntime

# get starting index for steady-state
it_min = ntime//3+1
it_max = ntime-1
# get number of time points in steady-state interval
it_interval = ntime - it_min

time_steady = time[it_min:it_max]
ntime_steady = time_steady.size

time_interval = time[it_max-1]-time[it_min]

def timeavg(ft):

    favg = simps(ft[it_min:it_max],x=time_steady) \
        / time_interval
    return favg
