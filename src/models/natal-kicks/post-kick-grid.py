
from pathlib import Path
import sys
import yaml

import poskiorb


def get_config_file(fname: str = "") -> dict:
    """Load configuration file with several options needed

    Parameters
    ----------
    fname : `string`
        Filename with config info

    Returns
    -------
    params : `dictionary`
        Info for the natal kicks run
    """

    configpath = Path(fname)
    if not configpath.is_file(): raise FileNotFoundError(f"`{fname}` not found.")

    try:
        with open(fname) as f:
            params = yaml.load(f, Loader=yaml.Loader)
    except Exception as e:
        raise e

    return params


def match_id(id: str="", fname: str="") -> float:
    """Search file for id and return value associated to it

    Parameters
    ----------
    id : string
        String to match

    Returns
    -------
    value: float
        Value that matches id
    """

    with open(fname, "r") as f:
        for k, line in enumerate(f):
            line = line.strip()
            if line.startswith(id):
                value = float(line.split(' ')[-1])
                return value


def get_collapse_data(star_fname: str, binary_fname: str) -> dict:
    """Retrieve information of collapsing binary, from mesabin2dco evolution

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
    """

    # check for existance of data files
    starpath = Path(star_fname)
    if not starpath.is_file(): raise FileNotFoundError(f"{star_fname} not found")

    binarypath = Path(binary_fname)
    if not binarypath.is_file(): raise FileNotFoundError(f"{binary_fname} not found")

    StarIdMaps = {
        "m1": "mass_pre_cc",
        "m1_core_mass": "c_core_mass_pre_cc",
        "m1_remnant_mass": "remnant_mass",
        "m1_fallback_fraction": "fallback_fraction"
    }
    BinaryIdMaps = {
        "m2": "companion_mass",
        "P": "period_pre_cc"
    }

    # get data pre-collapse
    CollapseData = dict()
    for key, value in StarIdMaps.items():
        CollapseData[key] = match_id(value, starpath)

    for key, value in BinaryIdMaps.items():
        CollapseData[key] = match_id(value, binarypath)

    return CollapseData


if __name__ == '__main__':

    # name of config file
    configname = sys.argv[1]

    # get configuration data
    ConfigData = get_config_file(fname=configname)

    # grab data at core collapse
    Data = get_collapse_data(ConfigData["Config"]["star-data-at-core-collapse"],
                             ConfigData["Config"]["binary-data-at-core-collapse"])


    # setup for natal kicks search
    binary = poskiorb.binary.BinarySystem(**Data)

    # reduce kick magnitudes by fallback if wanted
    if ConfigData["Config"]["reduce-kick-magnitude-with-fallback"]:
        f = lambda x: (1 - binary.m1_fallback_fraction) * x
    else:
        f = lambda x: x

    binary.set_natal_kick_distribution(n_trials=ConfigData["Config"]["number-of-kicks"],
                                       distribution_id=ConfigData["Config"]["kick-magnitude-distribution"],
                                       kick_scaling=f)

    # compute natal kicks & new orbital values
    binary.get_natal_kick_distribution()
    binary.get_orbital_distribution(verbose=True)
    binary.get_post_kick_grid(xnum=ConfigData["Config"]["xnum"],
                              ynum=ConfigData["Config"]["ynum"],
                              xquantiles=[ConfigData["Config"]["quantiles"]["xmin"], ConfigData["Config"]["quantiles"]["xmax"]],
                              
                              yquantiles=[ConfigData["Config"]["quantiles"]["ymin"], ConfigData["Config"]["quantiles"]["ymax"]],
                              min_prob=ConfigData["Config"]["min-prob"],
                              use_unbounded_for_norm=False, verbose=True)

    # show grid
    if ConfigData["Config"]["show-orbital-grid"]:
        binary.show_post_kick_with_grid(xattr="P", yattr="e", s=0.01)

    # store info on grid study
    if ConfigData["Config"]["save-kicks-info"]:
        binary.save_complete_grid(kick_fname=ConfigData["Config"]["kicks-fname"],
                                  orbit_fname=ConfigData["Config"]["orbits-after-kick-fname"])
        binary.save_target_grid(fname=ConfigData["Config"]["target-orbits-after-kick-fname"])
