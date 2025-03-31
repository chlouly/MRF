import numpy as np
from globals import *
import matplotlib.pyplot as plt

default_flips = [90, 120, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
#default_flips = [90, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]

# 117 mG * 1 ms = 180 degree flip

def fse_pulsetrain(pw, TR, flips=default_flips, fill=True):
    cur = 0
    out = np.zeros((ntime, 2))
    
    # This will be one RF pulse
    block_len = int(np.ceil(TR / dt))
    block = np.zeros(block_len)
    p_len = int(np.ceil(pw / dt))

    scale = pi / (180 * gam * p_len * dt)

    # The first pulse is diferent than the rest
    out[int(p_len / 2) : int(3 * p_len / 2) , x] = flips[0] * scale
    cur = block_len

    # This is our pulse shape for one pulse. We are using hanning
    for flip in flips[1:len(flips)]:
        block[0:p_len] = flip * scale
        if (cur + block_len) >= ntime:
            out[cur : ntime, y] = block
            return out
        else:
            out[cur : cur + block_len, y] = block
            cur = cur + block_len
    
    # If we get here, the loop didn't end early due to the ps filling up
    # If we need to fill then we do it here:
    while fill & ((cur + block_len) <= ntime):
        out[cur : cur + block_len, y] = block
        cur = cur + block_len

    #FIX ^^^

    return out