_DEFAULT_TAG = 'NULL'
__version__ = '0.0a0'


import datetime
import h5py
import numpy as np
import obspy
import os


def get_date_range(starttime, endtime):
    '''
    Get a generator of dates in a certain range.
    '''
    t_start = datetime.datetime(starttime.year, starttime.month, starttime.day)
    t_end = datetime.datetime(endtime.year, endtime.month, endtime.day)
    time_delta = t_end - t_start
    return ((starttime + datetime.timedelta(days=nday) for nday in range(0, time_delta.days+1)))


def get_first_samp_time(starttime, sampling_rate):
    '''
    Get the time of the first sample of the day.
    '''
    ts0 = obspy.UTCDateTime(f'{starttime.year}{starttime.julday:03d}')
    first_samp_time = starttime - int((starttime - ts0) * sampling_rate) / sampling_rate
    return (first_samp_time)


def get_sample_idx(time, starttime, sampling_rate, right=False):
    '''
    Get the index of a sample at a given time, relevant to starttime.
    '''
    if right is True:
        idx = int(np.ceil((time - starttime) * sampling_rate))
    else:
        idx = int(np.floor((time - starttime) * sampling_rate))
    return (idx)


class H5Seis(object):
    
    def __init__(self, filename, mode='a', compression='gzip', compression_opts=4):
        self._filename = filename
        self._create_dataset_kwargs = dict(compression=compression, compression_opts=compression_opts)
        self._h5 = h5py.File(filename, mode)
        for group in ('Waveforms',):
            if group not in self._h5.keys():
                self._h5.create_group(group)


    def add_waveforms(self, obj, tag=_DEFAULT_TAG):
        print(obj)
        if isinstance(obj, str) and os.path.isfile(os.path.abspath(obj)):
            st = obspy.read(obj)
        elif isinstance(obj, obspy.Stream):
            st = obj
        elif isinstance(obj, obspy.Trace):
            st = obspy.Stream(obj)
        else:
            raise(TypeError)
        self._add(st)


    def _add(self, st, tag=_DEFAULT_TAG):
#         if tag not in self._h5['/Waveforms'].keys():
#             self._h5.create_group(f'/Waveforms/{tag}')
# TODO: The Stream should be cleaned up here
        for tr in st:
            stats         = tr.stats
            network       = stats.network
            station       = stats.station
            location      = '__' if stats.location == '' else stats.location
            channel       = stats.channel
            starttime     = stats.starttime
            endtime       = stats.endtime
            delta         = stats.delta
            sampling_rate = stats.sampling_rate
            npts          = stats.npts
            key = f'/Waveforms/{tag}/{network}/{station}/{location}/{channel}'\
                  f'/{starttime.year}/{starttime.julday:03d}/{starttime}__{endtime}'
            if key not in self._h5:
                self._h5.create_dataset(
                    key, (npts,), 
                    dtype=tr.data.dtype,
                    **self._create_dataset_kwargs
                )
                self._h5[key].attrs['sampling_rate'] = sampling_rate
            self._h5[key][:] = tr.data


    def get_waveforms(self, network, station, location, channel, starttime, endtime, tag=_DEFAULT_TAG):
        key_base = f'/Waveforms/{tag}/{network}/{station}/{location}/{channel}'
        st = obspy.Stream()
        for date in get_date_range(starttime, endtime):
            day_key = f'{key_base}/{date.year}/{date.julday:03d}'
            if day_key not in self._h5:
                continue
            for item in self._h5[day_key]:
                ts, te = [obspy.UTCDateTime(s) for s in item.split('__')]
                if ts < endtime and te > starttime:
                    key = f'{day_key}/{item}'
                    sampling_rate = self._h5[key].attrs['sampling_rate']
                    idx_start = 0 if ts > starttime else get_sample_idx(starttime, ts, sampling_rate, right=True)
                    idx_end = None if te < endtime else get_sample_idx(endtime, ts, sampling_rate)
                    if idx_start == idx_end:
                        continue
                    data = self._h5[key][idx_start: idx_end]
                    tr = obspy.Trace(data=data)
                    tr.stats.starttime = ts + idx_start / sampling_rate
                    tr.stats.sampling_rate = sampling_rate
                    tr.stats.network       = network
                    tr.stats.station       = station
                    tr.stats.location      = '' if location == '__' else location
                    tr.stats.channel       = channel
                    st.append(tr)
        return (st)


    def close(self):
        self._h5.close()
