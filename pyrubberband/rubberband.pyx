import numpy as np

from libcpp cimport bool
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef enum Option:
    OptionProcessOffline       = 0x00000000
    OptionProcessRealTime      = 0x00000001

    OptionStretchElastic       = 0x00000000
    OptionStretchPrecise       = 0x00000010

    OptionTransientsCrisp      = 0x00000000
    OptionTransientsMixed      = 0x00000100
    OptionTransientsSmooth     = 0x00000200

    OptionDetectorCompound     = 0x00000000
    OptionDetectorPercussive   = 0x00000400
    OptionDetectorSoft         = 0x00000800

    OptionPhaseLaminar         = 0x00000000
    OptionPhaseIndependent     = 0x00002000

    OptionThreadingAuto        = 0x00000000
    OptionThreadingNever       = 0x00010000
    OptionThreadingAlways      = 0x00020000

    OptionWindowStandard       = 0x00000000
    OptionWindowShort          = 0x00100000
    OptionWindowLong           = 0x00200000

    OptionSmoothingOff         = 0x00000000
    OptionSmoothingOn          = 0x00800000

    OptionFormantShifted       = 0x00000000
    OptionFormantPreserved     = 0x01000000

    OptionPitchHighSpeed       = 0x00000000
    OptionPitchHighQuality     = 0x02000000
    OptionPitchHighConsistency = 0x04000000

    OptionChannelsApart        = 0x00000000
    OptionChannelsTogether     = 0x10000000

ctypedef int Options

cdef enum PresetOption:
    DefaultOptions             = 0x00000000
    PercussiveOptions          = 0x00102000

cdef extern from "RubberBandStretcher.h" namespace "RubberBand":
    cdef cppclass RubberBandStretcher:
        RubberBandStretcher(size_t, size_t, Options, double, double)
        void reset()
        void setTimeRatio(double ratio)
        void setPitchScale(double scale)
        double getTimeRatio() const
        double getPitchScale() const
        size_t getLatency() const
        void setTransientsOption(Options options)
        void setDetectorOption(Options options)
        void setPhaseOption(Options options)
        void setFormantOption(Options options)
        void setPitchOption(Options options)
        void setExpectedInputDuration(size_t samples)
        void setMaxProcessSize(size_t samples)
        size_t getSamplesRequired() const
        void setKeyFrameMap(const map[size_t, size_t] &)
        void study(const float *const *input, size_t samples, bool final)
        void process(const float *const *input, size_t samples, bool final)
        int available() const
        size_t retrieve(float *const *output, size_t samples) const
        float getFrequencyCutoff(int n) const
        void setFrequencyCutoff(int n, float f)
        size_t getInputIncrement() const
        vector[int] getOutputIncrements() const
        vector[float] getPhaseResetCurve() const
        vector[int] getExactTimePoints() const
        size_t getChannelCount() const
        void calculateStretch()
        void setDebugLevel(int level)

cdef extern from "main.h":
    string rubberband_process(vector[string] args, size_t frames, int samplerate, int channels, string input)


cdef __rubberband(y, sr, rbargs, verbose=False):

    ndim = y.ndim
    if ndim == 1:
        y = y[:, np.newaxis]

    frames, channels = y.shape

    rbargs = ['rubberband'] + rbargs

    if not verbose:
        rbargs.append('-q')

    rbargs = [str(arg).encode('utf-8') for arg in rbargs]

    dtype = y.dtype
    y = y.astype(np.float32)

    output = rubberband_process(rbargs, frames, sr, channels, y.tobytes('C'))
    result = np.frombuffer(output, dtype=np.float32)

    size = len(result)
    result = np.reshape(result, (size // channels, channels))

    if ndim == 1:
        result = np.squeeze(result, axis=1)

    return result.astype(dtype)


def time_stretch(y, sr, rate, rbargs=None, verbose=False):
    '''Apply a time stretch of `rate` to an audio time series.

    This uses the `tempo` form for rubberband, so the
    higher the rate, the faster the playback.


    Parameters
    ----------
    y : np.ndarray [shape=(n,) or (n, c)]
        Audio time series, either single or multichannel

    sr : int > 0
        Sampling rate of `y`

    rate : float > 0
        Desired playback rate.

    rbargs : list[str]
        Additional parameters for rubberband

        See `rubberband -h` for details.

    Returns
    -------
    y_stretch : np.ndarray
        Time-stretched audio

    Raises
    ------
    ValueError
        if `rate <= 0`
    '''

    if rate <= 0:
        raise ValueError('rate must be strictly positive')

    if rate == 1.0:
        return y

    if rbargs is None:
        rbargs = []

    rbargs = ['--tempo', rate] + rbargs

    return __rubberband(y, sr, rbargs, verbose=verbose)


def pitch_shift(y, sr, n_steps, rbargs=None, verbose=False):
    '''Apply a pitch shift to an audio time series.

    Parameters
    ----------
    y : np.ndarray [shape=(n,) or (n, c)]
        Audio time series, either single or multichannel

    sr : int > 0
        Sampling rate of `y`

    n_steps : float
        Shift by `n_steps` semitones.

    rbargs : list[str]
        Additional keyword parameters for rubberband

        See `rubberband -h` for details.

    Returns
    -------
    y_shift : np.ndarray
        Pitch-shifted audio
    '''

    if n_steps == 0:
        return y

    if rbargs is None:
        rbargs = []

    rbargs = ['--pitch', n_steps] + rbargs

    return __rubberband(y, sr, rbargs, verbose=verbose)
