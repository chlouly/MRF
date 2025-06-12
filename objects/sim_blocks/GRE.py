##########################################################################
#   This file contains all functions and classes related to simulating   #
#   the effects of FSE pulses.                                           #
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


class GRE(SimObj):
    """
    GRE(SimObj)

    This child class of SimObj represents a block in which we play
    a gradient echo.

    Inside of our ROI (the brain), the gradient echo readout is on resonance,
    therefor we set the semisolid absorption coefficient to 1.0 to represent
    perfect absorption.
    """
    absorption = 1.0  


    def  __init__(self, PW, ETL, delay, ESP, dt, sample_times=np.array([])):
        # This is the time of the block
        T = delay + (ETL * ESP)

        if (PW > ESP):
            raise ValueError("Error: Pulse width must be <= Echo Spacing ")             
        elif (T <= 0):
            raise ValueError("Error: The provided timing parameters created a block with 0 or negative time.")
        
        super().__init__(T, PW, ETL, delay, ESP, dt, sample_times=np.array([]))


    def set_rf(self, params):
        """
        This method overrides SimObj's definition of set_gradients(). Here we use
        the fse_pulsetrain() function to create a gradient echo RF pulsetrain.
        """

        # Get an array of the x and y components of the RF pulsetrain
        rf = gre_pulsetrain(self.PW, self.ESP, self.ETL, self.delay, self.T, self.dt, params.flip)

        # We will sample after the first echo for the center of k space
        RO_samples = np.array([self.delay + self.PW + self.dt])

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

    def set_flip(self, flip, phase=0):
        """
        This method lets you change the flip angle and RF phase of this block. It re-calls
        the method that sets the RF series.
        """

        # Reset the Effective B-Field
        self.B = np.zeros((self.ntime, 3))

        self.set_rf(flip, phase)
        self.set_gradients()    # Currently does nothing
    

# 117 mG * 1 ms = 180 degree flip
def gre_pulsetrain(PW, ESP, ETL, delay, T, dt, flip, phase=0):
    ntime = np.int32(np.ceil(T / dt))

    p_len = np.int32(np.ceil(PW / dt))
    d_len = np.int32(np.ceil(delay / dt))
    esp_len = np.int32(np.ceil(ESP / dt))

    # The block that represents one echo
    block = np.zeros((esp_len, 3))

    scale = flip * np.pi / (180 * gam * p_len * dt)
    block[0:p_len, 0] = scale * np.cos(phase)   # x component
    block[0:p_len, 1] = scale * np.sin(phase)   # y component

    # Repeat each echo the appropriate number of times
    out = np.tile(block, (ETL, 1))

    # adding delay at the beginning (if any)
    if d_len > 0:
        out = np.append(np.zeros((d_len, 3)), out, axis=0)

    # Now we pad the end (or check if weve gon over time)
    cur_out_len = np.shape(out)[0]

    if cur_out_len < ntime:
        # We have more time in the block so we play nothing
        out = np.append(out, np.zeros((ntime - cur_out_len, 3)), axis=0)
    elif cur_out_len > ntime:
        # Pulse is too long
        raise ValueError("Error: The given parameters generated an RF pulsetrain that does not fit in the aloughted block time.")
    
    return out
