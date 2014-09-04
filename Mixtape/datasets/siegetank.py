"""Siegetank integration

This module provides integration with siegetank.
"""
# Author: John Chodera <john.chodera@choderalab.org>
# Contributors:
# Copyright (c) 2014, Stanford University and the Authors
# All rights reserved.
#
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
from __future__ import print_function, absolute_import, division

from glob import glob
from io import BytesIO
from os import makedirs
from os.path import exists
from os.path import join
from zipfile import ZipFile
try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen

import mdtraj as md
from mixtape.datasets.base import Bunch
from mixtape.datasets.base import get_data_home

def strip_waters(trajectory):
    """Strip waters from a trajectory and return the stripped trajectory.

    Parameters
    ----------
    trajectory : Trajectory
        The trajectory from which waters are to be stripped.

    Returns
    -------
    stripped_trajectory : Trajectory
        A form of the input trajectory with waters stripped away.

    """

    atoms_to_keep = top.index[top.chainID == 0].values # TODO: assumption is that everything else is solvent?
    stripped_trajectory = trajectory.atom_slice(atoms_to_keep)
    return stripped_trajectory

def fetch_siegetank_target(target_token, login_token, data_home=None, download_if_missing=True, sync_seeds=False):
    """Retrieve data from SiegeTank SCV.

    Parameters
    ----------
    target_token : string
        Target token used to identify target from which streams are downloaded.

    login_token : string
        Token used to identify user for permission to log in.

    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all mixtape data is stored in '~/mixtape_data' subfolders.

    download_if_missing: optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    sync_seeds : optional, False by default
        If True, will sync initial seeds as well.

    Notes
    -----
    The siegetank package must be installed.
    """

    import siegetank

    # Log in using user token.
    siegetank.login(login_token)

    # Get target identifier.
    target = siegetank.load_target(target_token)

    # Construct description.
    description = "SiegeTank target %s\n" % target_token
    description += target.options['description']

    # Get list of streams.
    streams = target.streams()

    # Create directory to store mixtape data if it doesn't exist.
    data_home = get_data_home(data_home=data_home)
    if not exists(data_home):
        makedirs(data_home)

    # Create directory to store this stream if it doesn't exist.
    target_dir = join(data_home, target_token)
    if not exists(target_dir):
        print('Creating directory to store siegetank target %s : %s' % (target_token, target_dir))
        makedirs(target_dir)

    # Attempt to load all streams.
    trajectories = list()
    for stream in streams:
        # Get stream info.
        stream = siegetank.load_stream(stream_token)
        # Construct local data directory name for this data to go.
        data_folder = os.path.join(target_dir, stream.id)
        if not exists(data_folder):
            # Folder does not exist; sync data.
            makedirs(data_folder)
            stream.sync(data_folder, sync_seeds=sync_seeds)
            # Strip waters.
            stripped_trajectory = strip_waters(trajectory)
            # Write stripped trajectory over existing file.
            stripped_trajectory.save(filename)
        else:
            # Folder exists.
            # Load topology file.
            top = md.load(join(data_folder, 'reference.pdb'))
            trajectory = md.load(data_folder, 'stream.dcd'))
            # Load trajectory stripped of water.
            trajectory = md.load(join(data_folder, 'stream.dcd'))
            # Fetch more data from SCV if there are additional frames.
            if (stream.frames > trajectory.frames):
                # Fetch more data.
                # TODO
                pass
        # Append trajectory.
        trajectories.append(trajectory)


    # TODO: Strip out water.

    top = md.load(join(data_dir, 'ala2.pdb'))
    trajectories = []
    for fn in glob(join(data_dir, 'trajectory*.dcd')):
        trajectories.append(md.load(fn, top=top))

    # Return trajectories.
    return Bunch(trajectories=trajectories, DESCR=description)

fetch_alanine_dipeptide.__doc__ += __doc__

