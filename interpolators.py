import numpy as np
from scipy.interpolate import interp1d


# 'nearest' 
def interp_nearest(x_in, y_in, x_out, axis=0):
    f_ = interp1d(x_in, y_in, 
                  axis=axis, kind='nearest', fill_value='extrapolate')
    return f_(x_out)


# 'linear' 
def interp_linear(x_in, y_in, x_out, axis=0):
    f_ = interp1d(x_in, y_in, 
                  axis=axis, kind='linear', fill_value='extrapolate')
    return f_(x_out)


# 'cubic'
def interp_cubic(x_in, y_in, x_out, axis=0):
    f_ = interp1d(x_in, y_in, 
                  axis=axis, kind='cubic', fill_value='extrapolate')
    return f_(x_out)


# 'fftshift' (adapted from code provided by Ian Czekala)
def interp_fftshift(x_in, y_in, dx):
    # force even number of input channels
    Nch = y_in.shape[0]
    if (Nch % 2 == 1):
        if len(y_in.shape) > 1:
            y_add = np.atleast_3d(y_in[-1,:,:])
            y_in = np.append(y_in, np.transpose(y_add, (2, 0, 1)), axis=0)
        else:
            y_add = y_in[-1]
            y_in = np.append(y_in, y_add)
        x_in = np.append(x_in, x_in[-1] + (x_in[-1] - x_in[-2]))
        Nch += 1
        input_odd = True
    else:
        input_odd = False
    assert Nch % 2 == 0, "Only even-numbered input arrays."

    # input channel spacing
    dx_in = np.mean(np.diff(x_in))

    # Cross-correlation function using FFT
    rho_packed = np.fft.fft(np.fft.fftshift(y_in, axes=0), axis=0)
    fs_packed = np.fft.fftfreq(Nch, d=dx_in)

    # Apply the phase-shift
    phase = np.exp(-2.0j * np.pi * fs_packed * dx)
    if len(y_in.shape) > 1:
        phase = np.tile(phase, (y_in.shape[2], y_in.shape[1], 1)).T
    rho_packed_shifted = rho_packed * phase

    # inverse FFT back to frequency-space
    y_out = np.fft.fftshift(np.fft.ifft(rho_packed_shifted, axis=0), axes=0)

    # clip if padded
    if input_odd: y_out = y_out[:-1]

    return y_out



# =====

def interp_fftshift1d(xs, vals, nu_shift):
    nch = len(vals)
    if (nch % 2 == 1):
        vals = vals[:-1]
        xs = xs[:-1]
        nch -= 1
    assert nch % 2 == 0, "Only even numbered arrays for now."
    dchan = xs[1] - xs[0]

    # use fft to access cross-correlation function
    rho_packed = np.fft.fft(np.fft.fftshift(vals))
    fs_packed = np.fft.fftfreq(nch, d=dchan)

    # mulitply rho_packed by phase shift 
    # if shifting by positive dnu, then 
    phase = np.exp(-2.0j * np.pi * fs_packed * nu_shift)
    
    # transform back using ifft
    rho_packed_shifted = rho_packed * phase
    return np.fft.fftshift(np.fft.ifft(rho_packed_shifted))

def interp_fftshift_ORIG(xs, vals, nu_shift):
    """
    Shift the function by applying a phase shift in the frequency domain.
    This works for multi-dimensional vals with phase shift applied at
    the first dimension.

    Args:
        nu_shift: how much to shift by, in units of Hz
    """
    N = vals.shape[0]
    print(vals.shape)
    if (N % 2 == 1):
        vadd = np.atleast_3d(vals[-1,:,:])
        vals = np.append(vals, np.transpose(vadd, (2, 0, 1)), axis=0)
        xs = np.append(xs, xs[-1] + (xs[-1] - xs[-2]))
        N += 1
        input_odd = True
    else:
        input_odd = False
    assert N % 2 == 0, "Only even numbered arrays for now."
    dchan = xs[1] - xs[0]

    # use fft to access cross-correlation function
    rho_packed = np.fft.fft(np.fft.fftshift(vals, axes=0), axis=0)
    fs_packed = np.fft.fftfreq(N, d=dchan)

    # mulitply rho_packed by phase shift 
    # if shifting by positive dnu, then 
    phase = np.exp(-2.0j * np.pi * fs_packed * nu_shift)
    phase2 = np.tile(phase, (vals.shape[2], vals.shape[1], 1)).T

    # transform back using ifft
    rho_packed_shifted = rho_packed * phase2
    out = np.fft.fftshift(np.fft.ifft(rho_packed_shifted, axis=0), axes=0)

    if input_odd == True:
        out = out[:-1,:]
    return out
