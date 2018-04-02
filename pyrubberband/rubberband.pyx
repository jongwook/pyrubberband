from cpython.bytes cimport PyBytes_AsString, PyBytes_FromStringAndSize

import argparse
import numpy as np
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

__all__ = ['time_stretch', 'pitch_shift']

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


cdef __parse_args(rbargs):
    parser = argparse.ArgumentParser(description='rubberband arguments')
    parser.add_argument('--time', '-t', dest='time', type=float, help='Stretch to X times original duration')
    parser.add_argument('--tempo', '-T', dest='tempo', type=float, help='Change tempo by multiple X (same as --time 1/X)')
    parser.add_argument('--duration', '-D', dest='duration', type=float, help='Stretch or squash to make output file X seconds long')
    parser.add_argument('--pitch', '-p', dest='pitch', type=float, help='Raise pitch by X semitones')
    parser.add_argument('--frequency', '-f', dest='frequency', type=float, help='Change frequency by multiple X')

    parser.add_argument('--crisp', '-c', dest='crisp', type=int, help='Crispness (N = 0,1,2,3,4,5,6); default 5')
    parser.add_argument('--formant', '-F', dest='formant', action='store_true', help='Enable formant preservation when pitch shifting')

    parser.add_argument('--loose', '-L', dest='loose', action='store_true', help='Relax timing in hope of better transient preservation')
    parser.add_argument('--precise', '-P', dest='precise', action='store_true', help='Ignored: The opposite of -L, this is default from 1.6')
    parser.add_argument('--realtime', '-R', dest='realtime', action='store_true', help='Select realtime mode (implies --no-threads)')
    parser.add_argument('--no-threads', dest='no_threads', action='store_true', help='No extra threads regardless of CPU and channel count')
    parser.add_argument('--threads', dest='threads', action='store_true', help='Assume multi-CPU even if only one CPU is identified')
    parser.add_argument('--no-transients', dest='no_transients', action='store_true', help='Disable phase resynchronisation at transients')
    parser.add_argument('--bl-transients', dest='bl_transients', action='store_true', help='Band-limit phase resync to extreme frequencies')
    parser.add_argument('--no-lamination', dest='no_lamination', action='store_true', help='Disable phase lamination')
    parser.add_argument('--window-long', dest='window_long', action='store_true', help='Use longer processing window (actual size may vary)')
    parser.add_argument('--window-short', dest='window_short', action='store_true', help='Use shorter processing window')
    parser.add_argument('--smoothing', dest='smoothing', action='store_true', help='Apply window presum and time-domain smoothing')
    parser.add_argument('--detector-perc', dest='detector_perc', action='store_true', help='Use percussive transient detector (as in pre-1.5)')
    parser.add_argument('--detector-soft', dest='detector_soft', action='store_true', help='Use soft transient detector')
    parser.add_argument('--pitch-hq', dest='pitch_hq', action='store_true', help='In RT mode, use a slower, higher quality pitch shift')
    parser.add_argument('--centre-focus', dest='centre_focus', action='store_true', help='Preserve focus of centre material in stereo (at a cost in width and individual channel quality)')

    parser.add_argument('--debug', '-d', dest='debug', type=int, help='Select debug level (N = 0,1,2,3); default 0, full 3 (N.B. debug level 3 includes audible ticks in output)')
    parser.add_argument('--quiet', '-q', dest='quiet', action='store_true', help='Suppress progress output')

    return parser.parse_args(rbargs)


cdef __rubberband_params(args, sr, frames):
    cdef int options = 0
    cdef float time_ratio = 1.0
    cdef float pitch_scale = 1.0

    if args.time:
        time_ratio *= args.time
    if args.tempo:
        time_ratio /= args.tempo
    if args.duration:
        time_ratio = args.duration * sr / frames
    if args.pitch:
        pitch_scale *= 2.0 ** (args.pitch / 12.0)
    if args.frequency:
        pitch_scale *= args.frequency

    if args.crisp:
        if args.crisp == 0:
            args.no_transients = args.no_lamination = args.window_long = True
        if args.crisp == 1:
            args.detector_soft = args.no_lamination = args.window_long = True
        if args.crisp == 2:
            args.no_transients = args.no_lamination = True
        if args.crisp == 3:
            args.no_transients = True
        if args.crisp == 4:
            args.bl_transients = True
        if args.crisp == 6:
            args.no_lamination = args.window_short = True
    if args.formant:
        options |= RubberBandOptionFormantPreserved

    if not args.loose:
        options |= RubberBandOptionStretchPrecise
    if args.realtime:
        options |= RubberBandOptionProcessRealTime
    if args.no_threads:
        options |= RubberBandOptionThreadingNever
    if args.threads:
        options |= RubberBandOptionThreadingAlways
    if args.no_transients:
        options |= RubberBandOptionTransientsSmooth
    if args.bl_transients:
        options |= RubberBandOptionTransientsMixed
    if not args.no_lamination:
        options |= RubberBandOptionPhaseIndependent
    if args.window_long:
        options |= RubberBandOptionWindowLong
    if args.window_short:
        options |= RubberBandOptionWindowShort
    if args.smoothing:
        options |= RubberBandOptionSmoothingOn
    if args.detector_perc:
        options |= RubberBandOptionDetectorPercussive
    if args.detector_soft:
        options |= RubberBandOptionDetectorSoft
    if args.pitch_hq:
        options |= RubberBandOptionPitchHighQuality
    if args.centre_focus:
        options |= RubberBandOptionChannelsTogether

    return options, time_ratio, pitch_scale


cdef __rubberband(y, sr, rbargs, verbose=False):
    ndim = y.ndim
    if ndim == 1:
        y = y[:, np.newaxis]

    frames, channels = y.shape

    if not verbose:
        rbargs.append('-q')

    args = __parse_args([str(arg) for arg in rbargs])
    options, time_ratio, pitch_scale = __rubberband_params(args, sr, frames)

    dtype = y.dtype
    y = y.astype(np.float32)
    input_channels = [y[:, c].tobytes('C') for c in range(channels)]  # to retain ref count

    cdef float** input_data = <float **>malloc(channels * sizeof(float *))
    for c in range(channels):
        input_data[c] = <float *>PyBytes_AsString(input_channels[c])

    cdef RubberBandState rubberband = rubberband_new(sr, channels, options, time_ratio, pitch_scale)

    rubberband_set_expected_input_duration(rubberband, frames)
    rubberband_study(rubberband, input_data, frames, True)
    rubberband_process(rubberband, input_data, frames, True)

    cdef int available = rubberband_available(rubberband)
    output = np.zeros((channels, available), dtype=np.float32)
    cdef float [:, ::1] output_view = output
    cdef float **output_data = <float **>malloc(channels * sizeof(float *))

    cdef size_t ch = 0
    while ch < channels:
        output_data[ch] = <float *>&output_view[ch, 0]
        ch += 1

    cdef int remaining = available
    cdef int retrieved = 0
    while remaining > 0:
        retrieved = rubberband_retrieve(rubberband, output_data, remaining)
        remaining -= retrieved
        if remaining > 0:
            for c in range(channels):
                output_data[c] += retrieved

    rubberband_delete(rubberband)

    free(output_data)
    free(input_data)

    return output.astype(dtype).transpose()


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
