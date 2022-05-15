import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy import signal
from scipy.ndimage import convolve1d

characters = list(range(65, 91)) + list(range(48, 58))


def cut_epochs(t, data, timestamps):
    """
    cutting raw data into epochs
    :param t: tuple (start, end, samplerate)
    :param data: ndarray (n_channels, n_times)
    :param timestamps: list of timestamps
    :return: None or ndarray (n_epochs, n_channels, n_times)
    """
    assert data.ndim == 2
    timestamps = np.array(timestamps)
    start = timestamps + int(t[0] * t[2])
    end = timestamps + int(t[1] * t[2])
    epochs = np.stack([data[:, s:e] for s, e in zip(start, end)], axis=0)
    return epochs


def sort_epochs(epochs, event):
    """
    sorting epoch data according to event order
    :param epochs: 3d epochs (n_epochs, n_channels, time_seq)
    :param event: 2d event array (n_trials, 12)
    :return:
        sorted_epochs: ndarray, with the same shape of "epochs"
    """
    assert epochs.ndim == 3
    rep_dim = event.shape[1]
    indices = np.argsort(event, axis=-1).flatten()
    for i in range(0, indices.shape[0], rep_dim):
        indices[i:i + rep_dim] += i
    sorted_epochs = epochs[indices]  # deep copy
    return sorted_epochs


def timewindow(t_base, t_feat, epochs):
    """
    Cut specified time window from epoch data
    :param t_base: tuple, original time window (start, end)
    :param t_feat: tuple, target time window
    :param epochs: 3d array
    :return:
    """
    assert t_feat[0] >= t_base[0] and t_feat[1] <= t_base[1], print(t_feat, t_base)
    start, end = t_base
    unit = epochs.shape[-1] / (end - start)
    ind_s, ind_e = int((t_feat[0] - start) * unit), int((t_feat[1] - start) * unit)
    return epochs[..., ind_s:ind_e]


def common_average(data, ch_names=None):
    """

    :param data: (n_channels, n_times)
    :param ch_names: list
    :return:
    """
    average_data = data.mean(axis=0)
    filtered = data - average_data
    if ch_names is not None:
        ch_names = [ch + '_ave' for ch in ch_names]
        return filtered, ch_names
    else:
        return filtered


def bipolar(data, ch_names=None):
    """

    :param data: (n_channels, n_times)
    :param ch_names: list
    :return:
    """
    weight = np.array([1, -1], dtype=np.float64)[:, None]
    filtered = signal.convolve(data, weight, mode='valid')
    if ch_names is not None:
        ch_names_new = []
        for i in range(len(ch_names) - 1):
            ch_names_new.append(ch_names[i] + '_' + ch_names[i + 1])
        return filtered, ch_names_new
    else:
        return filtered


def laplacian(data, ch_names=None):
    """
    preprocessing
    :param data: data array 2D (n_channels, n_times)
                2d for seeg channels
    :param ch_names: list of channel names
    :return:
    """
    assert data.ndim == 2
    # 1D filter
    weight = np.array([1, -2, 1], dtype=np.float64)
    filtered = convolve1d(data, weight, axis=0, mode='reflect')
    if ch_names is not None:
        ch_names = [ch + '_lap' for ch in ch_names]
        return filtered, ch_names
    else:
        return filtered


def get_label(stim_string, n_rep):
    keyboard = {name: i for i, name in enumerate(characters)}
    assert isinstance(n_rep, int) and n_rep > 0

    label = np.zeros((len(stim_string), 12), dtype=np.int32)
    index = [keyboard[ord(char)] for char in stim_string]

    row_order = list(map(lambda x: x // 6, index))
    col_order = list(map(lambda x: x % 6 + 6, index))

    x = np.arange(len(stim_string))
    label[x, row_order] = 1
    label[x, col_order] = 1

    # tile
    label = np.tile(label[:, None, :], (1, n_rep, 1)).reshape((-1, 12))

    return label


########################################################################################################################
# utils for visualization
def plot_average(t, epoch, label, ch_names=None, cls_name=('nontarget', 'target'), with_error=True, **kwargs):
    """
    Draw average response
    :param t: tuple (start, end, samplerate)
    :param epoch: (epochs, channels, time)
    :param ch_names: list names of selected channels
    :param cls_name: list of classes name
    :param label: label (epochs, )
    :param with_error: whether or not to draw error shadow
    :param fig_size: figure size.
    :param colors: optional, iterable object of colors for each classes, str or matplotlib support colors
    :param fig: plt.figure object
    :param ax: plt.axes object
    :return: plt.figure object
    """
    start, end, fs = t
    t = np.linspace(start, end, epoch.shape[-1])
    unique_class = np.unique(label)
    cls_list = list(map(lambda i: epoch[label == i], unique_class))

    if ch_names is None:
        ch_names = [str(i) for i in range(epoch.shape[1])]

    if with_error:
        e_list = list(map(lambda x: np.std(x, axis=0) / np.sqrt(x.shape[0]), cls_list))

    colors = kwargs.pop('colors', None)
    alpha = kwargs.pop('alpha', 0.1)
    axes = kwargs.pop('axes', None)

    ave_cls_list = list(map(lambda x: np.mean(x, axis=0), cls_list))

    if axes is None:
        if epoch.shape[1] == 1:
            fig, axes = plt.subplots(1, 1, **kwargs)
            axes = np.array([axes])
        else:
            fig, axes = plt.subplots(int(math.ceil(len(ch_names) / 2)), 2, **kwargs)
    else:
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        fig = plt.gcf()

    for i, ax in enumerate(axes.flat):
        if i < len(ch_names):
            for j, (ave, cls) in enumerate(zip(ave_cls_list, cls_name)):
                if colors is not None:
                    ax.plot(t, ave[i], linestyle='-', label=cls, color=colors[j])
                else:
                    ax.plot(t, ave[i], linestyle='-', label=cls)

            if with_error:
                for j, (ave, e_ave, cls) in enumerate(zip(ave_cls_list, e_list, cls_name)):
                    if colors is not None:
                        ax.fill_between(t, ave[i] - e_ave[i], ave[i] + e_ave[i], alpha=alpha, color=colors[j])
                    else:
                        ax.fill_between(t, ave[i] - e_ave[i], ave[i] + e_ave[i], alpha=alpha)
            if i == 0:
                ax.legend()
            ax.set_title(ch_names[i], fontsize=14)
            ax.axvline(0, color='k')
            ax.tick_params('both', labelsize='large')
            ax.set_xticks([t[0], 0, t[-1]])
            ax.set_xticklabels(['%.2f' % t[0], '0', '%.2f' % t[-1]])
            ax.locator_params(axis='y', nbins=6)

    return fig


def draw_bar(data, err, bold, sig_channel, actual_channel, channel_type, show_legend=False):
    c = {'out': np.ones(3),
         'sec': '#96d7fa',
         'best': '#46a2f9'}
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 3))  # 60 mm
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    bars = plt.bar(sig_channel, data, yerr=err, width=0.8, color=[c[t] for t in channel_type],
                   edgecolor='black',
                   capsize=2)
    plt.tick_params('both', labelsize='large')
    plt.ylabel('Sig (SEEG HG)', fontsize=12)
    plt.xlabel('Channel', fontsize=12)
    #     plt.text(np.median(channel), ax.get_ylim()[1], subject, fontsize=14, horizontalalignment='center')
    plt.xticks(actual_channel)
    plt.xlim(0.5, 12.5)

    # bold data axis
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    bold, = ax2.plot(list(range(1, len(bold) + 1)), bold, linestyle='--', color='#E34239', alpha=1)
    ax2.set_ylabel('Sig (BOLD)', fontsize=12)
    ax2.tick_params('both', labelsize='large')

    # set legend
    if show_legend:
        out_index = channel_type.index('out')
        in_mt_index = channel_type.index('sec')
        best_index = channel_type.index('best')
        plt.legend((bars[best_index], bars[in_mt_index], bars[out_index], bold),
                   ('BEST', 'SEC', 'OUT', 'BOLD'), fontsize=12, frameon=False)
    plt.tight_layout()
    return fig
