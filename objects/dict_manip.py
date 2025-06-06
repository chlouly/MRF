import numpy as np
import h5py

# Globals - We use these to predefine the names of the fields stored in
# the dictionary file. (mostly to avoid bugs)
dict_name = "dictionary"        # Name of the field of dictionary entries in the file
idx_name = "cur_index"          # The name of the current parameter index (stored in case of crash)
T1f_name = "T1_f_vals"          # Array of T1 values that were simulated
T2f_name = "T2_f_vals"          # Array of T1 values that were simulated
T1s_name = "T1_s_vals"          # Array of T1 values that were simulated
alpha_name = "alpha_vals"       # Array of labeling efficiency values that were simulated
F_name = "F_vals"               # Array of perfusion rates that were simulated
ks_name = "ks_vals"             # Array of magnetization transfer rates that were simulated
kf_name = "kf_vals"             # Array of magnetization transfer rates that were simulated
CBV_name = "CBV_vals"           # Array of CBV values that were simulated
BAT_name = "BAT_vals"           # Array of BAT values that were simulated


def init_dict(name, params, num_samples):
    """
    This function creates and initializes a file in the HDF5 format that will store parameter values,
    dictionary entries, and the last-stored parameter indices (in case of crash).
    """
    if name == None:
        # If this happens, we will not be saving data, so we do nothing.
        return
    
    # Now we will initialize the file containing the dictionary.
    #try:
    d = h5py.File(name, "w")
    # Load the physiological parameter values into the file
    d.create_dataset(T1f_name, np.shape(params.T1_f_vals), data=params.T1_f_vals)
    d.create_dataset(T2f_name, np.shape(params.T2_f_vals), data=params.T2_f_vals)
    d.create_dataset(T1s_name, np.shape(params.T1_s_vals), data=params.T1_s_vals)
    d.create_dataset(alpha_name, np.shape(params.alpha_vals), data=params.alpha_vals)
    d.create_dataset(F_name, np.shape(params.F_vals), data=params.F_vals)
    d.create_dataset(ks_name, np.shape(params.ks_vals), data=params.ks_vals)
    d.create_dataset(kf_name, np.shape(params.kf_vals), data=params.kf_vals)
    d.create_dataset(CBV_name, np.shape(params.CBV_vals), data=params.CBV_vals)
    d.create_dataset(BAT_name, np.shape(params.BAT_vals), data=params.BAT_vals)

    # Set asside space for the last-stored parameter indices
    d.create_dataset(idx_name, np.size(params.get_shape()))
    
    # Set asside space for the actual dictionary
    #dict.create_dataset(dict_name, params.get_shape() + (num_samples,), compression='gzip')
    d.create_dataset(dict_name, params.get_shape() + (num_samples,))
    #dict.visititems(lambda name, obj: print(name))
    d.close()
    # except Exception as e:
    #     raise e


def store_entry(name, param_idx, entry):
    """
    This function stores an entry in the dictionary and updates the last simulated index.
    """

    with h5py.File(name, "r+") as dict:
        dict[dict_name][param_idx + (slice(None),)] = entry
        dict[idx_name][:] = list(param_idx)
