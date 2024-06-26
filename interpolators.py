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
