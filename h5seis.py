_DEFAULT_TAG     = 'NULL'
_LABEL_SEPARATOR = ';;'
__version__      = '0.0a0'


import datetime
import h5py
import itertools
import numpy as np
import obspy
import os
import re


def _filter(regex, iterable):
    '''
    Filter an iterable based on regular expression matches.
    '''
    return (filter(lambda key: re.match(regex, key), iterable))

def _cascade(keys, iterable):
    '''
    Cascade a set of regex filters for each level of /Waveforms/{tag}
    '''
    if len(keys) == 0:
        return (iterable)
    return ([_cascade(keys[1:], iterable[key]) for key in _filter(keys[0], iterable)])

def _flatten(iterable, n=1):
    '''
    Flatten the complex-shaped return value of _cascade()
    '''
    for i in range(n):
        iterable = itertools.chain.from_iterable(iterable)
    return (iterable)


def _get_date_range(starttime, endtime):
    '''
    Get a generator of dates in a certain range.
    '''
    t_start = datetime.datetime(starttime.year, starttime.month, starttime.day)
    t_end = datetime.datetime(endtime.year, endtime.month, endtime.day)
    time_delta = t_end - t_start
    return ((starttime + datetime.timedelta(days=nday) for nday in range(0, time_delta.days+1)))


#def _get_first_samp_time(starttime, sampling_rate):
#    '''
#    Get the time of the first sample of the day.
#    '''
#    ts0 = obspy.UTCDateTime(f'{starttime.year}{starttime.julday:03d}')
#    first_samp_time = starttime - int((starttime - ts0) * sampling_rate) / sampling_rate
#    return (first_samp_time)


def _get_sample_idx(time, starttime, sampling_rate, right=False):
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


    def __enter__(self):
        return (self)


    def __exit__(self, exc_type, exc_value, traceback):
        self._h5.close()


    @property
    def waveform_tags(self):
        return (sorted(list(self._h5['/Waveforms'])))


    def _add(self, st, tag=_DEFAULT_TAG, labels=None):
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
                ds = self._h5.create_dataset(
                    key, (npts,), 
                    dtype=tr.data.dtype,
                    **self._create_dataset_kwargs
                )
                ds.attrs['starttime']     = starttime.timestamp
                ds.attrs['sampling_rate'] = sampling_rate
                ds.attrs['network']       = network
                ds.attrs['station']       = station
                ds.attrs['location']      = location
                ds.attrs['channel']       = channel
                if labels is not None:
                    ds.attrs['labels'] = _LABEL_SEPARATOR.join(labels)
            else:
                ds = self._h5[key]
            ds[:] = tr.data

    
    def _get_trace(self, key, starttime=None, endtime=None):
        '''
        Extract a trace from dataset at {key} between {starttime} and {endtime}
        '''
        ts, te = [obspy.UTCDateTime(s) for s in key.split('/')[-1].split('__')]
        if starttime is not None and ts >= endtime:
            return (None)
        if endtime is not None and te <= starttime:
            return (None)
        ds = self._h5[key]
        sampling_rate = ds.attrs['sampling_rate']
        if starttime is None or ts > starttime:
            idx_start = 0
        else:
            idx_start = _get_sample_idx(starttime, ts, sampling_rate, right=True)
        if endtime is None or te < endtime:
            idx_end = None
        else:
            idx_end = _get_sample_idx(endtime, ts, sampling_rate)
        if idx_start == idx_end:
            return (None)
        data = ds[idx_start: idx_end]
        tr = obspy.Trace(data=data)
        tr.stats.starttime     = obspy.UTCDateTime(ds.attrs['starttime'] + idx_start / sampling_rate)
        tr.stats.sampling_rate = sampling_rate
        tr.stats.network       = ds.attrs['network']
        tr.stats.station       = ds.attrs['station']
        tr.stats.location      = '' if ds.attrs['location'] == '__' else ds.attrs['location']
        tr.stats.channel       = ds.attrs['channel']
        if 'labels' in ds.attrs:
            tr.stats.labels = ds.attrs['labels'].split(_LABEL_SEPARATOR)
        return (tr)
    
    
    def add_waveforms(self, obj, tag=_DEFAULT_TAG, labels=None):
        if isinstance(obj, str) and os.path.isfile(os.path.abspath(obj)):
            st = obspy.read(obj)
        elif isinstance(obj, obspy.Stream):
            st = obj
        elif isinstance(obj, obspy.Trace):
            st = obspy.Stream(obj)
        else:
            raise(TypeError)
        self._add(st, tag=tag, labels=labels)

        
    def close(self):
        self._h5.close()


    def get_waveforms(
        self, starttime, endtime,
        tag=_DEFAULT_TAG, network='.*', station='.*', location='.*', channel='.*'
    ):

        st = obspy.Stream()
        keys = (network, station, location, channel)
        groups = _cascade(keys, self._h5[f'/Waveforms/{tag}'])
        for group in _flatten(groups, n=3):
            for date in _get_date_range(starttime, endtime):
                day_key = f'{group.name}/{date.year}/{date.julday:03d}'
                if day_key not in self._h5:
                    continue
                for item in self._h5[day_key]:
                    key = f'{day_key}/{item}'
                    tr = self._get_trace(key, starttime=starttime, endtime=endtime)
                    if tr is not None:
                        st.append(tr)
        return (st)


    def get_waveforms_for_tag(
        self, tag,
        starttime=None, endtime=None, network='.*', station='.*', location='.*', channel='.*'
    ):
        st = obspy.Stream()
        keys = (network, station, location, channel)
        groups = _cascade(keys, self._h5[f'/Waveforms/{tag}'])
        for group in _flatten(groups, n=3):
            for year in group:
                for julday in group[year]:
                    for item in group[f'{year}/{julday}']:
                        key = f'{group.name}/{year}/{julday}/{item}'
                        tr = self._get_trace(key)
                        if tr is not None:
                            st.append(tr)
        return (st)