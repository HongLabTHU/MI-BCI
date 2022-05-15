import numpy as np
from mne.time_frequency import tfr_array_morlet
from scipy import signal


class FeatExtractor:
    def __init__(self, sfreq, band_erp=None, band_hg=None, method='fir', n=None):
        """
        Constructor of class featExtractor
        Pay attention: set erp_band to (1, ..) may result in really bad fir filter or extremely unstable iir filter
        So I recommend using lowpass filter to extract erp...
        Run testOfflineUtils to see if your modifications can pass the test.
        :param sfreq: int sample frequency
        :param band_erp: tuple f band of erp
        :param band_hg: tuple f band of hg
        :param method: 'iir' or 'fir'
        :param n: order of fir filter
        """
        assert method in ['iir', 'fir']
        assert band_erp is not None or band_hg is not None
        self.use_erp = band_erp is not None
        self.use_hg = band_hg is not None
        if band_erp is not None:
            band_erp = np.array(band_erp, dtype=np.float64)
            self._band_erp = band_erp
        if band_hg is not None:
            band_hg = np.array(band_hg, dtype=np.float64)
            self._band_hg = band_hg
        self.sfreq = sfreq

        if method == 'iir':
            if band_erp is not None:
                if band_erp.size == 1:
                    self.filter_erp = signal.butter(3, band_erp / sfreq * 2, btype='lowpass')
                else:
                    self.filter_erp = signal.butter(3, band_erp / sfreq * 2, btype='bandpass')
        else:
            if n is None:
                n = int(sfreq / 5)

            if band_erp is not None:
                if band_erp.size == 1:
                    b_erp = signal.firwin(n + 1, band_erp, fs=sfreq, pass_zero=True)
                else:
                    b_erp = signal.firwin(n + 1, band_erp, fs=sfreq, pass_zero=False)
                self.filter_erp = (b_erp, 1)

    def __call__(self, data, type):
        """
        extract features
        :param data: ndarray with the last axis "timesteps"
        :param type: "erp" or "hg"
        :return:
        """
        assert type in ['erp', 'hg']
        assert (type == 'erp' and self.use_erp) or (type == 'hg' and self.use_hg)
        if type == 'erp':
            erp = signal.filtfilt(*self.filter_erp, data, axis=-1)
            return erp
        elif type == 'hg':
            # morlet wavelet filter 10 Hz resolution
            freqs = np.arange(self._band_hg[0], self._band_hg[1], 10)
            n_cycles = freqs / 4
            power = tfr_array_morlet(data[None],
                                     sfreq=self.sfreq,
                                     freqs=freqs,
                                     n_cycles=n_cycles,
                                     output='avg_power',
                                     zero_mean=True)
            f_notch_ind = np.logical_not(np.isin(freqs, [50., 100., 150., 200.]))
            # remove power line noise, * f to normalize
            power = np.sum(power[:, f_notch_ind] * (freqs[f_notch_ind][None, :, None]), axis=1) / np.sum(
                freqs[f_notch_ind])
            hg_env = np.log10(power)
            hg_env = signal.detrend(hg_env, axis=-1, type='constant')
            return hg_env
