from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from objects import *

"""
This test is only to verify that the parameters are iterating correctly.
The correct order is:

    --FIRST--
    - CBV
    - ks
    - kf
    - T1_f
    - T2_f
    - T1_s
    - F
    - alpha
    - BAT
    --LAST--

Where the first ones are updated the most often and the last ones are
updated least often.
"""


if __name__ == "__main__":
    # We are only making sure that things are iterating correctly here
    # The input values don't really matter at all.
    tc = 1
    ta1 = np.array([6])
    ta2 = np.array([7, 8])
    ta3 = np.array([9, 10, 11])



    # Test a:
    # Checking to make sure single parameters still work as expected.
    print("\n\n-=-=- TEST 5a -=-=-")
    p1 = Params(tc, tc, tc, tc, tc, tc, tc, tc, tc, tc, tc, tc, tc, tc)
    p1_i = iter(p1)
    print("Makign sure things are initialized...")
    print("P1: CBV", p1.CBV)
    print("P1: BAT", p1.BAT)
    print("Now we will try looping over this item...")

    i = 0
    while True:
        print("\tLoop Index: ", i)
        i += 1
        print("\t", end="")
        p1_i.print_inds()
        try:
            next(p1_i)
        except StopIteration:
            break
    
    print("-=-=-=-=-=-=-=-=-=-")
    input("Press enter to continue with the Test 5...")



    # Test b:
    # Now we will throw in an array length 1 for ks, the behavior should remain the same
    print("\n\n-=-=- TEST 5b -=-=-")
    p2 = Params(tc, tc, tc, ta1, tc, tc, tc, tc, tc, tc, tc, tc, tc, tc)
    print("Makign sure things are initialized...")
    print("P1: CBV", p2.CBV)
    print("P1: BAT", p2.BAT)
    print("Now we will try looping over this item...")
    p2_i = iter(p2)

    i = 0
    while True:
        print("\tLoop Index: ", i)
        i += 1
        print("\t", end="")
        p2_i.print_inds()
        try:
            next(p2_i)
        except StopIteration:
            break

    print("-=-=-=-=-=-=-=-=-=-")
    input("Press enter to continue with the Test 5...")



    # Test c:
    # We will use an array of length 3 now for 2 parameters.
    print("\n\n-=-=- TEST 5c -=-=-")
    p3 = Params(tc, tc, ta3, tc, tc, tc, tc, ta3, tc, tc, tc, tc, tc, tc)
    print("Now we will try looping over this item...")

    p3_i = iter(p3)
    i = 0
    while True:
        print("\tLoop Index: ", i)
        i += 1
        print("\t", end="")
        p3_i.print_inds()
        try:
            next(p3_i)
        except StopIteration:
            break

    print("-=-=-=-=-=-=-=-=-=-")
    input("Press enter to continue with the Test 5...")



    # Test d:
    # We will use an array of length 3 now for 6 parameters.
    print("\n\n-=-=- TEST 5d -=-=-")
    p4 = Params(ta2, ta2, ta2, ta2, ta2, ta2, tc, tc, tc, tc, ta2, tc, tc, tc)
    print("Now we will try looping over this item...")
    p4_i = iter(p4)

    i = 0
    while True:
        print("\tLoop Index: ", i)
        i += 1
        print("\t", end="")
        p4_i.print_inds()
        try:
            next(p4_i)
        except StopIteration:
            break

    print("-=-=-=-=-=-=-=-=-=-")
    print("Test 5 complete!")

    