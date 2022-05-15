import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import f_classif

import utils
from dataset import Dataset
from model import FeatExtractor

# loading data
ch_names = [f'ch{i}' for i in range(2, 9)]
dataset = Dataset('S1', montage=['CH002_2K', 'CH003_2K', 'CH004_2K', 'CH005_2K', 'CH006_2K', 'CH007_2K', 'CH008_2K'])
print(dataset.data.shape)
print(len(dataset.timestamp))
print(dataset.events.shape)
print(ch_names)
fs = dataset.fs
down_sample_fs = 50.

target = utils.get_label('AHOV29', 10)
y = target.flatten()

# extract features
extractor = FeatExtractor(sfreq=dataset.fs,
                          band_erp=(1, 20),
                          band_hg=(60, 140),
                          method='fir')

# select channels in montage
data = dataset.data[dataset.montage_indices]

t_window = (-0.3, 0.9, fs)

# extract hg
hg = extractor(data, type='hg')
hg_epochs = utils.cut_epochs(t_window, hg, dataset.timestamp)
hg_epochs = utils.sort_epochs(hg_epochs, dataset.events)
down_ratio = int(fs / down_sample_fs)
hg_epochs_down_sampling = hg_epochs[..., ::down_ratio]

# extract erp
erp = extractor(data, type='erp')
erp_epochs = utils.cut_epochs(t_window, erp, dataset.timestamp)
erp_epochs = utils.sort_epochs(erp_epochs, dataset.events)
erp_epochs_down_sampling = erp_epochs[..., ::down_ratio]

# hg sig
hg_power_window = utils.timewindow(t_window[:2], (0, 0.5), hg_epochs_down_sampling)
hg_power_feature = hg_power_window.reshape((hg_power_window.shape[0], -1))
F, pval = f_classif(hg_power_feature, y)
sigval = -np.log10(pval)
sigval = sigval.reshape(hg_power_window.shape[1:])
fval = F.reshape(hg_power_window.shape[1:])

# draw figures
channel_type = ['out', 'best', 'best', 'best', 'sec', 'sec', 'out']
sig_channel = list(range(2, 9))
actual_channel = list(range(1, 9))
# bold sig values
bold = [0.04, 1.84, 8.63, 12.03, 8.94, 0.58, 0, 0]
subject = 'S1'
print(sigval.shape)
# top 25%
sig_top = np.sort(sigval, axis=1)[:, -5:]
utils.draw_bar(sig_top.mean(axis=1),
               sig_top.std(axis=1),
               bold, sig_channel,
               actual_channel,
               channel_type,
               show_legend=True)
# erp and hg average waveform
utils.plot_average(t_window, erp_epochs, y, ch_names, figsize=(10, 8))
utils.plot_average(t_window, hg_epochs, y, ch_names, figsize=(10, 8))
plt.show()
