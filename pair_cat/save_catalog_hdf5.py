import numpy as np
import warnings
import h5py
from utils import *

hh = 0.6774
c_mks = 3e8
msun_mks = 2e30
s_to_year = 3.17098e-8
year_to_s = 1.0 / s_to_year
lsun_ergs = 3.9e33
mdot_msun_yr = 1e10 / 980 / 1e6

METADATA = {
    # Header
    'SimulationName': None,
    'SimulationVersion': None,
    'NumberBinaries': None,
    'NumberBinariesDelay': None,
    'NumberAllBHs': None,
    'RedshiftsAllBHs': [None, 0.0, 40.0],
    'ModelType': None,
    'Version': None,
    'Date': None,
    'Contributor': None,
    'Email': None,
    'Principal': None,
    'Reference': None,
    'Website': None,

    # Model Parameters
    'HubbleConstant': ['km s^-1 Mpc^-1', 60.0, 80.0],
    'BoxSize': ['cMpc^3', 1.0e4, 1.0e9],
    'MinRedshift': [None, 0.0, 10.0],
    'MaxRedshift': [None, 0.0, 200.0],
    'MinBHSeedMass': ['M_sun', 1.0e2, 1.0e8],
    'StellarMassResolution': ['M_sun', 1.0e-1, 1.0e10],
    'DarkMatterMassResolution': ['M_sun', 1.0e-1, 1.0e12],
    'GasMassResolution': ['M_sun', 1.0e-1, 1.0e12],
    'SpatialResolution': ['kpc',0,1e2],
    'MinimumDarkMatterHaloMass': ['M_sun',1e4,1e15],
    'NoDelay': None,
    'Delay': None,
    'CommentsNoDelay': None,
    'CommentsDelay': None,
    'CommentsAllBHs': None,
}

# binaries
DATA_PAIRS = {
    "BlackHoles": {
        "PrimaryMass": ["M_sun", 1.0e2, 1.0e12],
        "SecondaryMass": ["M_sun", 1.0e2, 1.0e12],
        "Redshift": [None, 0.0, 40.0],
        "Separation": ["kpc", 1.0e-8, 1.0e5],
        "NumberDensity": ["Mpc^-3", 1e-30, 1e10],
    },
    "Galaxies": {
        "RemnantStellarMass": ["M_sun", 1.0e4, 1.0e14],
        "RemnantDarkMatterMass": ["M_sun", 1.0e4, 1.0e16],
        "RemnantRedshift": [None, 0.0, 40.0],
    },
}


# all BHs
DATA_ALLBHS = {
    "BlackHoles": {
        "Mass": ["M_sun", 1.0e2, 1.0e12],
        "Mdot": ["M_sun yr^-1", 0.0, 1.0e6],
        "RadEff": [None, 0.0, 10.0],
        "NumberDensity": ["Mpc^-3", 1e-30, 1e10],
    },
    "Galaxies": {
        "StellarMass": ["M_sun", 1.0e4, 1.0e14],
        "DarkMatterMass": ["M_sun", 1.0e4, 1.0e16],
    },
}




def write_hdf5(fname, metadata, singles, pairs, DEBUG):
    '''
    This function creates a hdf5 file containing the properties of the catalog. 
    The file has the following structure:
    - Header group

    - Pair group
        - NoDelay group
            - BlackHoles
                - PrimaryMass
                - SecondaryMass
                - Redshift
                - Separation
                - NumberDensity
            - Galaxies
                - RemnantStellarMass
                - RemnantDarkMatterMass
                - RemnantRedshift


    - All BHs group
        - Zxxx
            - BlackHoles
                - Mass
                - Mdot
                - RadEff
                - NumberDensity
            - Galaxies:
                - StellarMass
                - DarkMatterMass
    '''

    if DEBUG:
        print("\n###########################################")
        print("Writing hdf5 file:\n")

    # open file
    with h5py.File(fname, 'w') as h5out:
        # ---- write metadata "Header/" ----
        header = h5out.create_group("Header")
        u = header.create_group("units")
        for key, vals in metadata.items():
            if DEBUG:
                print(f"writing {key} = {vals}")
            header.attrs.create(key, data=vals)
            if isinstance(METADATA[key],list):
                if DEBUG:
                    print(f"writing units of {key} = {METADATA[key][0]}")
                u.attrs.create(key, data=str(METADATA[key][0]))           


        # ---- write pair data to "Pairs/" ----
        group = h5out.create_group("Pairs")
        for key, vals in pairs.items():
            if isinstance(vals, np.ndarray):
                if DEBUG:
                    print(f"writing {group.name}[{key}] = {type(vals)}, {vals.shape}, {vals.dtype}")
                tmp = group.create_dataset(key, data=vals)
            # else:
            #     tmp.attrs.create(key,data=value)

        # ---- write All BHs data to "AllBHs/" ----
        allbhs_group = h5out.create_group("AllBHs")
        for key, vals in singles.items():
            if isinstance(vals, np.ndarray):
                if DEBUG:
                    print(f"writing {group.name}[{key}] = {type(vals)}, {vals.shape}, {vals.dtype}")
                tmp = allbhs_group.create_dataset(key, data=vals)
            # else:
            #     tmp.attrs.create(key,data=value)
    h5out.close()
    if DEBUG:
        print(f"Done writing to {fname}")
        print("###########################################\n")

    return


