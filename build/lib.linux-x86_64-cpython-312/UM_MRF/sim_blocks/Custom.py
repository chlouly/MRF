##########################################################################
#   This file contains all functions and classes related to simulating   #
#   the effects of a custom pulse.                                       #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

from .SimObj import SimObj
from numpy import size, shape

class Custom(SimObj):
    """
    Custom(SimObj)

    This child class of SimObj is meant to be a blank slate that can take in any arbitrary
    effective B field and incorporate it in the simulation. At the moment, there is no input
    validation for the shape of B and s to amke sure they are arrays of the propper shape.
    """


    def __new__(cls, B, s, dt):
        if (shape(B)[1] != 3) | (len(shape(B)) > 2):
            raise ValueError("Error: B(t) must have 3 spatial dimensions (Must be of shape (n, 3))")
        if shape(s)[0] != shape(B)[0]:
            raise ValueError("Error: s(t) and B(t) must have the same number of timepoints")
        
        T = len(s) * dt

        return super().__new__(cls, T)
        

    def __init__(self, B, s, dt, dynamic_time=False):  
        """
        Creates an instance of the Custom class - Custom(SimObj)

        Input Arguments:
            B:          (n, 3) numpy array describing the effective B field in the rotating frame
                        (Including RF excitement and possibly gradients)
            s:          (n, 1) numpy array describing the arterial activation function.
            dt:         Timestep [ms]
        """
        ntime = len(s)
        T = ntime * dt

        super().__init__(T, 0, 0, 0, 0, dt, dynamic_time=dynamic_time)

        self.B = B
        self.s = s


    def set_rf(self, params):
        """
        This method overrides SimObj's definition of set_rf(). The RF pulses should already be
        incorporated into this object's B field array, therefor there is nothing to do.

        This method essentially does nothing.
        """
        pass


    def set_gradients(self):
        """
        This method overrides SimObj's definition of set_gradients(). The gradient fields should 
        already be incorporated into this object's B field array, therefor there is nothing to do.

        This method essentially does nothing.
        """
        pass

    