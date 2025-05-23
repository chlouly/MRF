##########################################################################
#   This file contains all functions and classes related to simulating   #
#   the effects of pCASL pulses.                                         #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

import numpy as np
from ..SimObj import SimObj

pi = np.pi
gambar = 42570                  # Gyromagnetic coefficient [kHz/T]
gam = gambar * 2 * 3.14159      # Gamma [kRad/sT]

x = 0
y = 1


class pCASL(SimObj):
    """
    pCASL(SimObj)

    This child class of SimObj represents a block in which we play
    a pCASL pulsetrain.

    Inside of our ROI (the brain), the fast spin echo readout is played off
    resonance, because of this, we expect the pulse to saturate the semisolid
    pool, and to not directly affect the tissue pool. We use a saturation term
    as a constant to set the rate of the exponential decay of the semisolid
    magnetization. This choice of constant is arbitrary.
    """
    absorption = 0.25
    saturation = 1000


    def __init__(self, T, dt, sample_times=np.array([]), control=0):
        """
        Creates an instance of the DeadAir class - DeadAir(SimObj)

        Input Arguments:
            params:     Instance of Params object
            T:          Block length [ms]
            dt:         Timestep [ms]
            control:    Bool for label vs control pulse. (default = 0)
        """
        self.control = control
        super().__init__(T, 0, 0, dt, sample_times)


    def set_rf(self):
        """
        This method overrides SimObj's definition of set_rf(). We don't directly
        simulate the effects of a pCASL pulse train, we merely simulate its expected
        effects on the two pools (nothing for the free pool, and saturation for the
        semisolid pool). Therefor we want no modification to this class' B field array 
        (it is initialized to 0s).

        This method essentially does nothing.
        """
        pass


    def set_gradients(self):
        """
        This method overrides SimObj's definition of set_gradients(). There should be no
        gradients playing since we are not directly simulating a pCASL pulse, therefor 
        we want no modification to this class' B field array (it is initialized to 0s).

        This method essentially does nothing.
        """
        pass


    def set_s_shape(self, time_queue, BAT):
        """
        This method overrides SimObj's definition of set_s_sig(). The arterial longitudinal
        magnetization is affected by a pCASL pulse, and for now we assume that it takes the
        the form of a rect function. Here, we create a start and end time for the rect function
        created by this pulse, we add those times to a queue in the form of a tuple, and we pass 
        the updated queue off to SimObj's definition of set_s_sig().

        See SimObj's definition of this method for details...
        """
        # We only see activation in the labeling case
        if not self.control:
            # The rect function starts BAT ms from the beginning of the pulse
            # The rect ends BAT ms after the the pulsetrain ends BAT + durration of pulse ms
            # from the beginning of the pulse. So we add:
            #       (Bolus Start [ms], Bolus End [ms])
            #    =  (BAT [ms], BAT + self.T [ms])
            time_queue += [(BAT, BAT + self.T)]
        
        # We call SimObj's definition of set_s_sig()
        return super().set_s_shape(time_queue, BAT)


def pcasl_rf_gen(flip, ntime, dt, pw, TR, d_psi, psi_0, control=False):
    # This will be one RF pulse
    block = np.zeros(np.int64(np.ceil(TR / dt)))
    p_len = int(np.ceil(pw / dt))

    # This is our pulse shape for one pulse. We are using hanning
    block[0:p_len] = np.hanning(p_len)

    H = np.sum(block[0:p_len]) * dt
    amp = flip * pi / (H * 180 * gam)

    block[0:p_len] = block[0:p_len] * amp

    b_len = np.size(block)

    # For control, we add phi to the phase 
    if control:
        d_psi = d_psi + pi
    
    # Create the final pulse sequence and return
    n_reps = np.int64(np.ceil(ntime / b_len))
    B_out = np.zeros((n_reps * b_len, 3))

    for i in range(n_reps):
        B_out[i * b_len : (i + 1) * b_len, x] = block * np.cos(psi_0 + i * d_psi)
        B_out[i * b_len : (i + 1) * b_len, y] = block * np.sin(psi_0 + i * d_psi)

    return B_out[0:ntime, :]


# def slice_grad_gen(Gmax, pw, Gave, TR, dt, ntime):
#     # This will be one Gradient pulse
#     block_len = int(np.ceil(TR / dt))
#     block = np.zeros(block_len)
#     p_len = int(np.ceil(pw / dt))

#     # Construct the ramp up to the max amplitude
#     main_ramp = np.arange(0, Gmax, max_slew * dt)
#     main_ramp_len = np.size(main_ramp)

#     # Construct the main gradient pulse for one repetition
#     block[0:main_ramp_len] = main_ramp
#     block[main_ramp_len:main_ramp_len + p_len] = Gmax
#     block[main_ramp_len + p_len:2 * main_ramp_len + p_len] = main_ramp[::-1]

#     # Construct the rephasing pulse
#     cur_g_sum = np.sum(block) * dt
#     reph_g_sum = Gave - cur_g_sum
#     reph_amp = np.sqrt(np.abs(reph_g_sum * max_slew * dt))
#     reph_ramp = np.arange(0, -reph_amp, -max_slew * dt)
#     reph_ramp_len = np.size(reph_ramp)

#     block[2 * main_ramp_len + p_len:2 * main_ramp_len + p_len + reph_ramp_len] = reph_ramp
#     block[2 * main_ramp_len + p_len + reph_ramp_len:2 * main_ramp_len + p_len + 2 * reph_ramp_len] = reph_ramp[::-1]

#     # Create the final pulse sequence and return
#     nreps = int(np.ceil(ntime / block_len))
#     return np.tile(block, (nreps))[0 : ntime]

