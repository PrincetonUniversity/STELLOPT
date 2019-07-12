# -*- coding: utf-8 -*-

# load currents in the coil system of W7-X from the archive, take mean over
# a given experiment ID timeslice and return the set of 7 numbers
# 2017/09/07, H. Thomsen: adaptions to new read_restdb
# (which now includes get_start_end)

from W7Xrest.read_restdb import read_restdb, get_start_end
import numpy as _np

# read all coil signals as one stream
#__stream_url = "http://archive-webapi.ipp-hgw.mpg.de/ArchiveDB/raw/W7X/CoDaStationDesc.84/DataReductionProcessDesc.30_DATASTREAM/_signal.json"

# fast time traces: 5kHz, +-1A
# but only absoulte values...
baseurl="http://archive-webapi.ipp-hgw.mpg.de/ArchiveDB/raw/W7X/CoDaStationDesc.84/"
#streamurls=["DataModuleDesc.20043_DATASTREAM/0/AAE10_Idx1-1",
#            "DataModuleDesc.20055_DATASTREAM/0/AAE29_Idx1-1",
#            "DataModuleDesc.20067_DATASTREAM/0/AAE38_Idx1-1",
#            "DataModuleDesc.20079_DATASTREAM/0/AAE47_Idx1-1",
#            "DataModuleDesc.20091_DATASTREAM/0/AAE56_Idx1-1",
#            "DataModuleDesc.20103_DATASTREAM/0/AAE14_Idx1-1",
#            "DataModuleDesc.20115_DATASTREAM/0/AAE23_Idx1-1"]
            
# fast time traces: 5kHz, +-1A
streamurls=["DataModuleDesc.21643_DATASTREAM/0/current_npc1 (from AAE10)",
            "DataModuleDesc.21643_DATASTREAM/1/current_npc2 (AAE29)",
            "DataModuleDesc.21643_DATASTREAM/2/current_npc3 (AAE38)",
            "DataModuleDesc.21643_DATASTREAM/3/current_npc4 (AAE47)",
            "DataModuleDesc.21643_DATASTREAM/4/current_npc5 (AAE56)",
            "DataModuleDesc.21643_DATASTREAM/5/current_pca (AAE14)",
            "DataModuleDesc.21643_DATASTREAM/6/current_pca (AAE23)"]


#%% fetch currents timetraces for given interval
def get_w7x_currents_timetrace(t_from, t_upto):
    import numpy as np
    # get start and end time of program    
    #[valid, times, current_data] = read_restdb(__stream_url+"?from="+_np.str(t_from)+"&upto="+_np.str(t_upto))
    all_currents=[]
    all_times=[]
    all_valid=True
    for i in range(len(streamurls)):
        #print(baseurl+streamurls[i])
        [valid,times,current_data]=read_restdb(baseurl+streamurls[i]+"/_signal.json?from="+_np.str(t_from)+"&upto="+_np.str(t_upto))
        if valid:
            all_times.append(times)
            all_currents.append(current_data)
        else:
            all_valid=False
    
    if not all_valid:
        print("error: could not get valid current data!")
        return [None, None]
    else:
        currents=np.mean(all_currents, axis=1)
        return [all_times[0], currents]
    return [None, None]

#%% fetch mean currents (and stddev if requested) for given experiment ID
def get_w7x_currents_byID(program_id):
    # get start and end time of program
    [t_from, t_upto] = get_start_end(program_id)
    import numpy as np

    currents=np.zeros([7,1])
    all_valid=True
    for i in range(len(streamurls)):
        #print(baseurl+streamurls[i])
        [valid,times,current_data]=read_restdb(baseurl+streamurls[i]+"/_signal.json?from="+_np.str(t_from)+"&upto="+_np.str(t_upto))
        if valid:
            #print(current_data)
            currents[i]=np.mean(current_data)
        else:
            all_valid=False
            
    if not all_valid:
        print("error: could not get valid current data!")
        return None
    else:
        return currents

#%% fetch mean currents (and stddev if requested) for given experiment ID
def get_w7x_currents(nanosecond):
    import numpy as np
    # fast measurement of main coil currents is done at a rate of 5kHz
    # => 200us between samples => using this we should always get only one sample, i.e. the nearest one
    t_from=nanosecond-100000    
    t_upto=nanosecond+100000-1
    
    currents=np.zeros([7,1])
    all_valid=True
    for i in range(len(streamurls)):
        #print(baseurl+streamurls[i])
        [valid,times,current_data]=read_restdb(baseurl+streamurls[i]+"/_signal.json?from="+_np.str(t_from)+"&upto="+_np.str(t_upto))
        if valid:
            print(current_data)
            currents[i]=current_data
        else:
            all_valid=False
            
    if not all_valid:
        print("error: could not get valid current data!")
        return None
    else:
        return [times, currents]

