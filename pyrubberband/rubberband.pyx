import numpy as np

cdef extern from "rubberband-c.h":
    cdef enum RubberBandOption:
        RubberBandOptionProcessOffline       = 0x00000000
        RubberBandOptionProcessRealTime      = 0x00000001

        RubberBandOptionStretchElastic       = 0x00000000
        RubberBandOptionStretchPrecise       = 0x00000010

        RubberBandOptionTransientsCrisp      = 0x00000000
        RubberBandOptionTransientsMixed      = 0x00000100
        RubberBandOptionTransientsSmooth     = 0x00000200

        RubberBandOptionDetectorCompound     = 0x00000000
        RubberBandOptionDetectorPercussive   = 0x00000400
        RubberBandOptionDetectorSoft         = 0x00000800

        RubberBandOptionPhaseLaminar         = 0x00000000
        RubberBandOptionPhaseIndependent     = 0x00002000

        RubberBandOptionThreadingAuto        = 0x00000000
        RubberBandOptionThreadingNever       = 0x00010000
        RubberBandOptionThreadingAlways      = 0x00020000

        RubberBandOptionWindowStandard       = 0x00000000
        RubberBandOptionWindowShort          = 0x00100000
        RubberBandOptionWindowLong           = 0x00200000

        RubberBandOptionSmoothingOff         = 0x00000000
        RubberBandOptionSmoothingOn          = 0x00800000

        RubberBandOptionFormantShifted       = 0x00000000
        RubberBandOptionFormantPreserved     = 0x01000000

        RubberBandOptionPitchHighSpeed       = 0x00000000
        RubberBandOptionPitchHighQuality     = 0x02000000
        RubberBandOptionPitchHighConsistency = 0x04000000

        RubberBandOptionChannelsApart        = 0x00000000
        RubberBandOptionChannelsTogether     = 0x10000000

    ctypedef int RubberBandOptions;
    ctypedef void* RubberBandState;

    RubberBandState rubberband_new(unsigned int sampleRate,
                                      unsigned int channels,
                                      RubberBandOptions options,
                                      double initialTimeRatio,
                                      double initialPitchScale);

    void rubberband_delete(RubberBandState);

    void rubberband_reset(RubberBandState);

    void rubberband_set_time_ratio(RubberBandState, double ratio);
    void rubberband_set_pitch_scale(RubberBandState, double scale);

    double rubberband_get_time_ratio(const RubberBandState);
    double rubberband_get_pitch_scale(const RubberBandState);

    unsigned int rubberband_get_latency(const RubberBandState);

    void rubberband_set_transients_option(RubberBandState, RubberBandOptions options);
    void rubberband_set_detector_option(RubberBandState, RubberBandOptions options);
    void rubberband_set_phase_option(RubberBandState, RubberBandOptions options);
    void rubberband_set_formant_option(RubberBandState, RubberBandOptions options);
    void rubberband_set_pitch_option(RubberBandState, RubberBandOptions options);

    void rubberband_set_expected_input_duration(RubberBandState, unsigned int samples);

    unsigned int rubberband_get_samples_required(const RubberBandState);

    void rubberband_set_max_process_size(RubberBandState, unsigned int samples);
    void rubberband_set_key_frame_map(RubberBandState, unsigned int keyframecount, unsigned int *, unsigned int *);

    void rubberband_study(RubberBandState, const float *const *input, unsigned int samples, int final);
    void rubberband_process(RubberBandState, const float *const *input, unsigned int samples, int final);

    int rubberband_available(const RubberBandState);
    unsigned int rubberband_retrieve(const RubberBandState, float *const *output, unsigned int samples);

    unsigned int rubberband_get_channel_count(const RubberBandState);

    void rubberband_calculate_stretch(RubberBandState);

    void rubberband_set_debug_level(RubberBandState, int level);
    void rubberband_set_default_debug_level(int level);


cdef __rubberband(y, sr, rbargs, verbose=False):

    ndim = y.ndim
    if ndim == 1:
        y = y[:, np.newaxis]

    frames, channels = y.shape

    rbargs = ['rubberband'] + rbargs

    if not verbose:
        rbargs.append('-q')

    rbargs = [str(arg).encode('utf-8') for arg in rbargs] # TODO argparse
    options = RubberBandOptionProcessOffline

    dtype = y.dtype
    y = y.astype(np.float32)

    rubberband = rubberband_new(sr, channels, options, 1.0, 1.0)
    rubberband_delete(rubberband)

    return 'NYI'


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
