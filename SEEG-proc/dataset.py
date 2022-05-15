import glob
import os

import numpy as np
import scipy.io as sio

from nex import Reader as NexReader


def _load_nex(data_dir):
    """
    nex file loader
    :param data_dir:
    :return:
        data: ndarray, shape (ch, timesteps)
        ch_names: list,  name of each channel
        timestamps: list, stimulation onset
    """
    files = glob.glob(os.path.join(data_dir, '*.nex'))
    assert len(files) == 1

    reader = NexReader(useNumpy=True)
    data = reader.ReadNexFile(files[0])

    var = data['Variables']
    ch_names = []
    trigger_ch = None
    con_data = []
    samplerate = None
    for i, ch in enumerate(var):
        if 'CH' in ch['Header']['Name']:
            ch_names.append(ch['Header']['Name'])
            con_data.append(ch['ContinuousValues'])
            samplerate = ch['Header']['SamplingRate']
        if 'digin' == ch['Header']['Name']:
            trigger_ch = i
    assert trigger_ch is not None
    timestamp = np.round(data['Variables'][trigger_ch]['Timestamps'] * samplerate).astype(np.int32).tolist()
    con_data = np.array(con_data)
    return samplerate, con_data, ch_names, timestamp


class Dataset:
    """
    for loading data and stimulation order.
    """
    data_format = {
        'nex': _load_nex,
    }

    def __init__(self, subject, montage=None):
        self.subject = subject
        self._subj_path = os.path.join(os.path.dirname(__file__), '../data/SEEG', subject)
        self.root_dir = os.path.join(os.path.dirname(__file__), '../data/SEEG', subject)

        self.montage = montage

        self.events = self.load_event()
        # load data and timestamps
        self.fs, dataarray, ch_names, timestamp = self._load_data()
        timestamp = Dataset.ts_check(timestamp, self.fs)
        self.data = dataarray
        # list to set
        self.ch_names = ch_names
        self.timestamp = timestamp
        self.montage_indices = self.get_channel_indices(self.montage, self.ch_names)

    def _load_data(self):
        """
        Read data according to file format
        :return:
            dataext: str, data file name

        """
        walk_path = self.root_dir
        loader = None
        for f in os.listdir(walk_path):
            _ext = f.split('.')[-1]
            try:
                loader = Dataset.data_format[_ext]
                break
            except KeyError:
                pass
        if loader is None:
            raise FileNotFoundError('No matching data format found')
        return loader(walk_path)

    def load_event(self):
        walk_path = self.root_dir
        file = glob.glob(os.path.join(walk_path, self.subject) + '*')
        assert len(file) == 1
        file = file[0]

        if file.endswith('.mat'):
            raw = sio.loadmat(file)
            order = raw['stim_order']
            order -= 1
            return order.reshape((-1, 12))
        else:
            with open(file) as f:
                stim_order = [[int(x) for x in line.split()] for line in f if len(line) > 1]
            return np.array(stim_order)

    @staticmethod
    def get_channel_indices(target_channels, channels_in_data):
        """
        Get corresponding index number for channels in target channels
        :param target_channels: list, target channel names
        :param channels_in_data: list, all channel names in data source.
        :return:
        """
        indices = []
        # build a dictionary for indexing
        channel_book = {name: i for i, name in enumerate(channels_in_data)}
        for ch in target_channels:
            try:
                indices.append(channel_book[ch])
            except ValueError as err:
                print(err)

        return indices

    @staticmethod
    def ts_check(ts, fs):
        # check the first time stamp. 
        # Sometimes there would be an additional wrong trigger in the beginning.
        while len(ts) % 12 and (not (fs * 0.1 <= ts[1] - ts[0] <= fs * 0.3)):
            del ts[0]
        return ts
