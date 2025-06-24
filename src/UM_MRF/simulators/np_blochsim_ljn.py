##########################################################################
#   This file contains all functions definitions used to run the Bloch   #
#   Flow simulation using Jaynes' Matrix formalism as outlined in the    #
#   2020 paper:                                                          #
#       "Numerical approximation to the general kinetic model            #
#       for ASL quantification"                                          #
#   written by:                                                          #
#       Nam G. Lee, Ahsan Javed, Terrence R. Jao, Krishna S. Nayak       #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

import numpy as np
import time as tm

# Constant
gambar = 42570                  # Gyromagnetic coefficient [kHz/T]
gam = gambar * 2 * 3.14159      # Gamma [kRad/sT]

def is_approx(a, b, rel_tol=1e-12, abs_tol=0.0):
    """
    Helper function to see if two floats a and b are approvimately equal
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def get_rot_mat(u, theta, B_mag, absorption):
    """
    This is a helpter function that creates a rotation matrix specifically
    to be used in the LJN algorithm. Elements (1:3, 1:3) of the matrix is 
    a rotation matrix of the following form:
        R(u, theta) rotates by theta radians about the unit vector u
    where:
        u = B / |B|
        theta = gamma |B| dt
    The (4, 4) element of the matrix holds the term that tips the semisolid
    longitudinal magnetization based on its absorption spectrum. This term is
    not very rigorously implemented, mostly because it is not too important to
    what we are trying to do.

    Parameters:
        u:          Unit vector about which we will rotate      (length 3)
        theta:      The number of radians that we will rotate by
        B_mag:      The L2 norm of the effective B field.
        absorption: The coefficient that governs the effect of the pulse
                    on the semisolid pool.

    Output:
        R:          A 4x4 matrix following the form described above.
    """
    # First, we're going to pre-calulate some values
    ux, uy, uz = u
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    one_minus_cos = 1 - cos_t

    # We form the final matrix and return it
    return np.array([
        # The first 3 rows correspond to the rotation matrix
        [cos_t + ux**2 * one_minus_cos,
        ux*uy*one_minus_cos - uz*sin_t,
        ux*uz*one_minus_cos + uy*sin_t,
        0.0],

        [uy*ux*one_minus_cos + uz*sin_t,
        cos_t + uy**2 * one_minus_cos,
        uy*uz*one_minus_cos - ux*sin_t,
        0.0],

        [uz*ux*one_minus_cos - uy*sin_t,
        uz*uy*one_minus_cos + ux*sin_t,
        cos_t + uz**2 * one_minus_cos,
        0.0],

        # This final row corresponds to the semisolid absorption
        [0.0, 0.0, 0.0, np.exp(-np.pi * (gam * B_mag)**2 * absorption)]
    ])


def np_ljn_setp(M, B, s, p, dt, ACE, absorption, s_sat):
    """
    Function that performs one step of the LJN Simulation

    Parameters:
        M:              The previous Magnetization vector M(t - dt)     (length 4)
        B:              The current effective B field vector B(t)       (length 3)
        s:              The current arterial magnetization s(t)
        p:              Instance of a Params object
        dt:             Timestep [ms]
        ACE:            The precomputed Matrix formed by the products of:
                            A(t), C(t) and E(t)
                        As described in the paper from which this method
                        was obtained.
        absorption:     The arbitrary constant that governs the absorption
                        of the semisolid pool.
        s_sat:          The arbitrary constant that only serves to saturate
                        the semisolid pool durring pCASL pulses (set to 0 otherwise)

    Output:
        out:            The simulated Magnetization vector at the current time M(t)
                        (length 4)
    """

    # The first step is to simulate the errects of the RF pulse and the Gradients
    # (if any) in the form of a rotation matrix around the 
    B_mag = np.linalg.norm(B)
    if is_approx(0.0, B_mag):
        # If the magnitude of our B field approximately 0, we don't have to
        # do anything.
        out = M
    else:
        # Otherwise, we find a unit vector about which we will rotate
        # and the angle theta
        u = B / B_mag
        theta = gam * B_mag * dt

        # Get calculate rotation matrix
        R = get_rot_mat(u, theta, B_mag, absorption)

        # Apply the rotation matrix
        out = np.matmul(R, M)


    # This term is here to artoficially saturate the semisolid pool.
    # If the saturation constant is not zero, the semisolid longitudinal
    # magnetization exponentially decays to 0, and quite quickly.
    # Otherwise, the semisolid pool is left alone
    #   
    # We only do this durring the pCASL pulses to force the expected behavior
    # without actually having to simulate a real pCASL pulse.
    if not is_approx(s_sat, 0):
        out[3] = out[3] * np.exp(-np.pi * s_sat)


    # Now we calculate the pseudo-steady-state term from the paper
    # We will call it D(t), even though it is not actually reffered to as D
    # in the original paper. This vector in the original paper is:
    #   (Lambda + Gamma + Xi)^-1 * D(t)
    D = np.zeros(4)
    D[2] = (1 + p.T1_s * p.ks) * (p.M0_f + s * p.T1f_app) + p.T1f_app * p.ks * p.M0_s
    D[3] = (p.T1_s * p.ks) * (p.M0_f + s  * p.T1f_app) + (1 + p.T1f_app * p.ks) * p.M0_s

    D = D / (1 + p.T1f_app * p.kf + p.T1_s * p.ks)


    # Now we perform the a simplification of the final step
    # to obtain our final output vector.
    # return M[t + 1] = ACE * M[t] + (I - ACE) * D
    #                 = ACE * (M[t] - D) + D
    return np.matmul(ACE, out - D) + D


def np_blochsim_ljn(B, s, p, dt, n_time, M_start, absorption, s_sat=0.0, crusher_inds=np.array([]), timer=False):
    """
    This function runs a simulation using the method outlined in the LJN paper.

    Parameters:
        B:              The (n_time, 3) array of effective B field values.
        s:              The (n_time, ) array of arterial magnetization values.
        dt:             Timestep [ms]
        n_time:         The number of timepoints that we are using in this simulation.
        M_start:        The initial Magnetization vector M[0, :]
        absorption:     The arbitrary constant that governs the absorption
                        of the semisolid pool.
        s_sat:          The arbitrary constant that only serves to saturate
                        the semisolid pool durring pCASL pulses (set to 0 otherwise)
        timer:          A boolean flag that indicates wether or not we wish to time
                        this simulation (default value is False)

    Output:
        M:              The (n_time, 4) array of simulated Magnetization vectors.
    """
    # Getting a start time in case we choose to time this simulation
    start_t = tm.time()


    # We precompute the product A(t) C(t) E(t) because it only depends
    # on the size of the timestep, not any parameters that change throughout
    # the simulation.
    # A(t) matrix
    A = np.identity(4, dtype=np.float64)
    A_expval = np.exp(-(1 + p.f) * p.ks * dt)
    A[2, 2] = (1 + p.f * A_expval) / (1 + p.f)
    A[3, 2] = (p.f - p.f * A_expval) / (1 + p.f)
    A[2, 3] = (1 - A_expval) / (1 + p.f)
    A[3, 3] = (p.f + A_expval) / (1 + p.f)

    # CE matrix (C times E)
    CE = np.identity(4, dtype=np.float64)
    CE[0, 0] = np.exp(-dt * p.R2f_app)
    CE[1, 1] = np.exp(-dt * p.R2f_app)
    CE[2, 2] = np.exp(-dt * p.R1f_app)
    CE[3, 3] = np.exp(-dt * p.R1s_app)

    # Not we multiply to get the final matrix
    ACE = np.matmul(A, CE)


    # Now we bedin the simulation
    # First, we allocate an array for the simulated Magnetization values
    # as well as setting it's initial value.
    M = np.zeros((n_time, 4))
    M[0, :] = M_start

    # Now we loop through all timepoints and run one step at each timepoint
    cur_crush_ind = 0
    for t in range(1, n_time):
        M[t, :] = np_ljn_setp(M[t - 1, :], B[t, :], s[t], p, dt, ACE, absorption, s_sat)
        if (np.size(crusher_inds) > 0) and (t == crusher_inds[cur_crush_ind]):
            # Crush it
            M[t, 0] = 0
            M[t, 1] = 0

            # Next index
            cur_crush_ind = (cur_crush_ind + 1) % np.size(crusher_inds)


    # We get the time at which the simulation ended in case we chose to time it
    end_t = tm.time()

    # If timer was true, we subtract the start from the end and display how long it took.
    if timer:
        print(f"Time: {end_t - start_t}")

    return M

