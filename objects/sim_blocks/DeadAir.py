##########################################################################
#   This file contains all functions and classes related to simulating   #
#   the effects of Dead Air (nothing playing on scanner).                #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

from ..SimObj import SimObj

# No Pulse Playing
class DeadAir(SimObj):
    """
    DeadAir(SimObj)

    This child class of SimObj represents a block where the scanner is not
    playing any RF excitations or gradients. 
    """


    def __init__(self, T, dt):
        """
        Creates an instance of the DeadAir class - DeadAir(SimObj)

        Input Arguments:
            params:     Instance of Params object
            dt:         Timestep [ms]
        """
        super().__init__(T, 0, 0, 0, 0, dt)


    def set_rf(self, params):
        """
        This method overrides SimObj's definition of set_rf(). There should be no
        RF pulses in a DeadAir block, therefor we want no modification to this
        class' B field array (it is initialized to 0s)

        This method essentially does nothing.
        """
        pass


    def set_gradients(self):
        """
        This method overrides SimObj's definition of set_gradients(). There should be no
        gradients playing in a DeadAir block, therefor we want no modification to this
        class' B field array (it is initialized to 0s)

        This method essentially does nothing.
        """
        pass
