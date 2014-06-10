"""Utility functions"""
# Author: Robert McGibbon <rmcgibbo@gmail.com>
# Contributors: Bharath Ramsundar <bharath.ramsundar@gmail.com>
# Copyright (c) 2014, Stanford University
# All rights reserved.

# Mixtape is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Mixtape. If not, see <http://www.gnu.org/licenses/>.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from __future__ import print_function, division, absolute_import

import os
import json
import hashlib
import numpy as np
from sklearn.utils import check_random_state
from sklearn.externals.joblib import load, dump
from numpy.linalg import norm
import mdtraj as md

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------


def sha1digest(filename):
    """SHA1 hex digest of a file
    """
    sha1 = hashlib.sha1()
    with open(filename, 'rb') as f:
        sha1.update(f.read())
    return sha1.hexdigest()


def verbosedump(value, fn, compress=1):
    """verbose wrapper around joblib.dump"""
    print('Saving "%s"... (%s)' % (fn, type(value)))
    dump(value, fn, compress=compress)


def verboseload(fn):
    """verbose wrapper around joblib.load"""
    print('loading "%s"...' % fn)
    return load(fn)


def verbosedump_or_skip(function, fn, compress=1, digest=None):
    if os.path.exists(fn):
        try:
            loaded = verboseload(fn)
            if digest is None or loaded['digest'] == digest:
                print('Skipping %s. Already exists' % fn)
                return
        except:
            pass

    value = function()
    verbosedump({'value': value,
                 'digest': digest},
                fn, compress=compress)


def iterobjects(fn):
    for line in open(fn, 'r'):
        if line.startswith('#'):
            continue
        try:
            yield json.loads(line)
        except ValueError:
            pass

def gridscores_to_dataframe(grid):
    """Convert the grid_scores_ from a DistributedGridSearchCV
    run into a pandas dataframe that is easier to plot and
    manipulate

    Parameters
    ----------
    grid_scores_ : list of dicts
    """
    import pandas
    items = []
    for elem in grid:
        item = {'test_mean' : elem['cv_validation_scores'].mean(),
                'test_std' : elem['cv_validation_scores'].std()}
        item.update(elem['parameters'])
        if 'cv_train_scores' in elem:
            item['train_mean'] = elem['cv_train_scores'].mean()
            item['train_std'] = elem['cv_train_scores'].std()
        items.append(item)

    return pandas.DataFrame(items)


def rmsd(X, Y, yi):
    # md.rsmd isn't picklable, so this is a little proxy
    return md.rmsd(X, Y, yi, precentered=True)


def categorical(pvals, size=None, random_state=None):
    """Return random integer from a categorical distribution

    Parameters
    ----------
    pvals : sequence of floats, length p
        Probabilities of each of the ``p`` different outcomes.  These
        should sum to 1.
    size : int or tuple of ints, optional
        Defines the shape of the returned array of random integers. If None
        (the default), returns a single float.
    random_state: RandomState or an int seed, optional
        A random number generator instance.
    """
    cumsum = np.cumsum(pvals)
    if size is None:
        size = (1,)
        axis = 0
    elif isinstance(size, tuple):
        size = size + (1,)
        axis = len(size) - 1
    else:
        raise TypeError('size must be an int or tuple of ints')

    random_state = check_random_state(random_state)
    return np.sum(cumsum < random_state.random_sample(size), axis=axis)


##########################################################################
# MSLDS Utils (experimental)
##########################################################################


def iter_vars(A, Q, N):
    """Utility function used to solve fixed point equation
       Q + A D A.T = D
       for D
     """
    V = np.eye(np.shape(A)[0])
    for i in range(N):
        V = Q + np.dot(A, np.dot(V, A.T))
    return V


##########################################################################
# END of MSLDS Utils (experimental)
##########################################################################

def map_drawn_samples(selected_pairs_by_state, trajectories):
    """Lookup trajectory frames using pairs of (trajectory, frame) indices.

    Parameters
    ----------
    selected_pairs_by_state : np.ndarray, dtype=int, shape=(n_states, n_samples, 2)
        selected_pairs_by_state[state, sample] gives the (trajectory, frame)
        index associated with a particular sample from that state.
    trajectories : list(md.Trajectory)
        The trajectories assocated with sequences,
        which will be used to extract coordinates of the state centers
        from the raw trajectory data

    Returns
    -------
    frames_by_state : mdtraj.Trajectory, optional
        If `trajectories` are provided, this output will be a list
        of trajectories such that frames_by_state[state] is a trajectory
        drawn from `state` of length `n_samples`
    
    Examples
    --------
    >>> selected_pairs_by_state = hmm.draw_samples(sequences, 3)
    >>> samples = map_drawn_samples(selected_pairs_by_state, trajectories)
    
    Notes
    -----
    YOU are responsible for ensuring that selected_pairs_by_state and 
    trajectories correspond to the same dataset!
    
    See Also
    --------
    utils.map_drawn_samples : Extract conformations from MD trajectories by index.
    ghmm.GaussianFusionHMM.draw_samples : Draw samples from GHMM    
    ghmm.GaussianFusionHMM.draw_centroids : Draw centroids from GHMM    
    """

    frames_by_state = []

    for state, pairs in enumerate(selected_pairs_by_state):
        frames = [trajectories[trj][frame] for trj, frame in pairs]
        state_trj = np.sum(frames)  # No idea why numpy is necessary, but it is
        frames_by_state.append(state_trj)
    
    return frames_by_state
