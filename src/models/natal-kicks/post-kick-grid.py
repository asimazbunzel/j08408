
from pathlib import Path
import sys

import poskiorb


def match_id(id: str='', fname: str=''):
    '''Search file for id and return value associated to it

    Parameters
    ----------
    id : string
        String to match

    Returns
    -------
    value: float
        Value that matches id
    '''

    with open(fname, 'r') as f:
        for k, line in enumerate(f):
            line = line.strip()
            if line.startswith(id):
                value = float(line.split(' ')[-1])
                return value


def get_collapse_data(star_fname: str, binary_fname: str):
    '''Retrieve information of collapsing binary, from mesabin2dco evolution
    
    Parameters
    ----------
    star_fname : string
        Filename with star data

    binary_fname : string
        Filename with binary data

    Returns
    -------
    CollapseData : dictionary
        Important data needed to compute a grid of post-kick orbital parameters
    '''

    # check for existance of data files
    starpath = Path(star_fname)
    if not starpath.is_file(): raise FileNotFoundError(f'{star_fname} does not exist')

    binarypath = Path(binary_fname)
    if not binarypath.is_file(): raise FileNotFoundError(f'{binary_fname} does not exist')

    StarIdMaps = {
        'm1': 'mass_pre_cc',
        'm1_core_mass': 'c_core_mass_pre_cc',
        'm1_remnant_mass': 'remnant_mass',
        'm1_fallback_fraction': 'fallback_fraction'
    }
    BinaryIdMaps = {
        'm2': 'companion_mass',
        'P': 'period_pre_cc'
    }

    # get data pre-collapse
    CollapseData = dict()
    for key, value in StarIdMaps.items():
        CollapseData[key] = match_id(value, starpath)

    for key, value in BinaryIdMaps.items():
        CollapseData[key] = match_id(value, binarypath)

    return CollapseData


if __name__ == '__main__':

    starname = sys.argv[1]
    binaryname = sys.argv[2]

    Data = get_collapse_data(starname, binaryname)

    binary = poskiorb.binary.BinarySystem(**Data)

    binary.set_natal_kick_distribution(n_trials=5000, distribution_id='Maxwell',
                                       kick_scaling=lambda x: (1-binary.m1_fallback_fraction)*x)

    binary.get_natal_kick_distribution()

    binary.get_orbital_distribution(verbose=True)

    binary.get_post_kick_grid(use_unbounded_for_norm=False, verbose=True)

    # binary.show_post_kick_with_grid(xattr='P', yattr='e', s=1)

    binary.save_complete_grid(kick_fname='kick-grid.data', orbit_fname='orbit.data')

    # binary.save_complete_grid(fname='kick-grid.data')
