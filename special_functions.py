#!/usr/bin/env python

# This is a python re-writing of Shanjie Zhang and Jianming Jin's code, as
# attached to their book "Computation of Special Functions." That code is
# explicitly licensed for use in other programs, provided acknowledgement that
# the original work is copyrighted:
#
# Copyright (c) 1996 John Wiley & Sons, Inc.

from math import log

def python_comelp(hk):

    pk = 1.0e0 - hk**2

    if hk == 1.0:
         ck = 1.0e300
         ce = 1.0
    else:
         ak = (((.01451196212*pk + .03742563713)*pk + .03590092383)*pk + .09666344259)*pk + 1.38629436112
         bk = (((.00441787012*pk + .03328355346)*pk + .06880248576)*pk + .12498593597)*pk + .5
         ck = ak - bk*log(pk)
         ae = (((.01736506451*pk + .04757383546)*pk + .0626060122)*pk + .44325141463)*pk + 1.0
         be = (((.00526449639*pk + .04069697526)*pk + .09200180037)*pk + .2499836831)*pk
         ce = ae - be*log(pk)

    return ck, ce
