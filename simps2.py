# Simps2.py: Modified simpson's method of solving an integral from
# scipy.integrate for integration in findBedLoad function
# This code is modified from SciPy.


# Copyright © 2001, 2002 Enthought, Inc.
# All rights reserved.
#
# Copyright © 2003-2019 SciPy Developers.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#     Neither the name of Enthought nor the names of the SciPy Developers may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


from __future__ import division, print_function, absolute_import

from scipy.special.orthogonal import p_roots
from scipy.special import gammaln
from numpy import sum, ones, add, diff, isinf, isscalar, \
     asarray, real, trapz, arange, empty, abs
# Added abs for calculating absolute value of differences
import numpy as np
import math
import warnings


def _basic_simps(y,start,stop,x,dx,axis):
    nd = len(y.shape)
    if start is None:
        start = 0
    step = 2
    all = (slice(None),)*nd
    slice0 = tupleset(all, axis, slice(start, stop, step))
    slice1 = tupleset(all, axis, slice(start+1, stop+1, step))
    slice2 = tupleset(all, axis, slice(start+2, stop+2, step))

    if x is None:  # Even spaced Simpson's rule.
        result = add.reduce(dx/3.0 * (y[slice0]+4*y[slice1]+y[slice2]),
                                    axis)
    else:
        # Account for possibly different spacings.
        #    Simpson's rule changes a bit.
        h = abs(diff(x,axis=axis))
        # Changed h to absolute of the difference account for moving along
        # the outside of a circle
        sl0 = tupleset(all, axis, slice(start, stop, step))
        sl1 = tupleset(all, axis, slice(start+1, stop+1, step))
        h0 = h[sl0]
        h1 = h[sl1]
        hsum = h0 + h1
        hprod = h0 * h1
        h0divh1 = h0 / h1
        result = add.reduce(hsum/6.0*(y[slice0]*(2-1.0/h0divh1) +
                                              y[slice1]*hsum*hsum/hprod +
                                              y[slice2]*(2-h0divh1)),axis)
    return result


def simps2(y, x=None, dx=1, axis=-1, even='avg'):
    y = asarray(y)
    nd = len(y.shape)
    N = y.shape[axis]
    last_dx = dx
    first_dx = dx
    returnshape = 0
    if x is not None:
        x = asarray(x)
        if len(x.shape) == 1:
            shapex = [1] * nd #ones(nd)
            shapex[axis] = x.shape[0]
            saveshape = x.shape
            returnshape = 1
            x = x.reshape(tuple(shapex))
        elif len(x.shape) != len(y.shape):
            raise ValueError("If given, shape of x must be 1-d or the "
                    "same as y.")
        if x.shape[axis] != N:
            raise ValueError("If given, length of x along axis must be the "
                    "same as y.")
    if N % 2 == 0:
        val = 0.0
        result = 0.0
        slice1 = (slice(None),)*nd
        slice2 = (slice(None),)*nd
        if even not in ['avg', 'last', 'first']:
            raise ValueError("Parameter 'even' must be 'avg', 'last', or 'first'.")
        # Compute using Simpson's rule on first intervals
        if even in ['avg', 'first']:
            slice1 = tupleset(slice1, axis, -1)
            slice2 = tupleset(slice2, axis, -2)
            if x is not None:
                last_dx = x[slice1] - x[slice2]
            val += 0.5*last_dx*(y[slice1]+y[slice2])
            result = _basic_simps(y,0,N-3,x,dx,axis)
        # Compute using Simpson's rule on last set of intervals
        if even in ['avg', 'last']:
            slice1 = tupleset(slice1, axis, 0)
            slice2 = tupleset(slice2, axis, 1)
            if x is not None:
                first_dx = x[tuple(slice2)] - x[tuple(slice1)]
            val += 0.5*first_dx*(y[slice2]+y[slice1])
            result += _basic_simps(y,1,N-2,x,dx,axis)
        if even == 'avg':
            val /= 2.0
            result /= 2.0
        result = result + val
    else:
        result = _basic_simps(y,0,N-2,x,dx,axis)
    if returnshape:
        x = x.reshape(saveshape)
    return result


def tupleset(t, i, value):
    l = list(t)
    l[i] = value
    return tuple(l)
