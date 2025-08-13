from .SimObj import SimObj
import numpy as np

# No Pulse Playing
class BIR8(SimObj):
    """
    BIR8(SimObj)

    This child class of SimObj represents a block where the scanner is not
    playing any RF excitations or gradients. 
    """

    eTE = 22
    crush_length = 750 


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

    
    def run_ljn(self, p, M_start=...):
        self.M = np.zeros((self.ntime, 4))

        # T2 decay
        self.M[:, 2] = M_start[2] * np.exp(-self.eTE / p.T2_f)

        # Crush
        self.M[:, 0:2] = 0.0

        # T1 decay curing crusher
        self.M[:, 2] *= -np.exp(-self.crush_length / p.T1_f)
        self.M[:, 2] += 1

        # I think this will saturate the semisolid pool
        self.M[:, 3] = 0.0

        return self.M


