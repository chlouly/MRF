##########################################################################
#   This file contains all functions and classes related to simulating   #
#   the effects of FSE pulses.                                           #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

import numpy as np
from globals import *
from ..SimObj import SimObj


class FSE(SimObj):
    """
    FSE(SimObj)

    This child class of SimObj represents a block in which we play
    a fast spin echo.

    Inside of our ROI (the brain), the fast spin echo readout is on resonance,
    therefor we set the semisolid absorption coefficient to 1.0 to represent
    perfect absorption.
    """
    absorption = 1.0    


    def set_rf(self):
        """
        This method overrides SimObj's definition of set_gradients(). Here we use
        the fse_pulsetrain() function to create a fast spin echo RF pulsetrain.
        """

        # Get an array of the x and y components of the RF pulsetrain
        rf = fse_pulsetrain(self.PW, self.TR, self.T, self.dt)

        # We will sample after the first 90y
        RO_samples = np.array([self.TR])

        # Append the ReadOut sample times to the sample times array
        self.sample_times = np.append(self.sample_times, RO_samples)

        # Call the parent class' definition of set_rf() to add the pulse to
        # the objects effective B field.
        super().set_rf(rf)


    def set_gradients(self):
        """
        This method overrides SimObj's definition of set_gradients(). At the moment
        we are not considering the effects of gradients, however this function was left
        in in case we decide to in the future.

        This method essentially does nothing.
        """
        pass
    



default_flips = [90, 120, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
#default_flips = [90, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
# 117 mG * 1 ms = 180 degree flip
def fse_pulsetrain(pw, TR, T, dt, flips=default_flips, fill=True):
    cur = 0
    ntime = np.int64(np.ceil(T / dt))
    out = np.zeros((ntime, 3))
    
    # This will be one RF pulse
    block_len = int(np.ceil(TR / dt))
    block = np.zeros(block_len)
    p_len = int(np.ceil(pw / dt))

    scale = pi / (180 * gam * p_len * dt)

    # The first pulse is diferent than the rest
    out[np.int64(p_len / 2) : np.int64(3 * p_len / 2) , x] = flips[0] * scale
    cur = block_len

    # This is our pulse shape for one pulse. We are using hanning
    for flip in flips[1:len(flips)]:
        block[0:p_len] = flip * scale
        if (cur + block_len) >= ntime:
            out[cur : ntime, y] = block[0 : ntime - cur]
            return out
        else:
            out[cur : cur + block_len, y] = block
            cur = cur + block_len
    
    # If we get here, the loop didn't end early due to the ps filling up
    # If we need to fill then we do it here:
    while fill & ((cur + block_len) < ntime):
        out[cur : cur + block_len, y] = block
        cur = cur + block_len

    return out
