"""LISA AstroWG - MBH Catalogs - Catalog generation script - v2.1

Catalog v2.1 is primarily meant to collect 
(i) information on binary MBHs for “no delay” cases that will be used to create a uniform “delay” case for all models, 
(ii) information on the full MBH population that will be used to build AGN luminosity functions and MBH-Mgal relations to compare to observations. 
A group can also provide a “delay” case that will be used for the paper, but this can also be provided at later time. 
There are some additional data, such as dark matter halo information, that can be optionally provided for v2.1 if desired.  

Links
-----
- Catalog v2.1 spreadsheet:
  https://docs.google.com/spreadsheets/d/1Zz-A5QcKHdPabpSmkKcs9ERzeAA0sqM_/edit?usp=sharing&ouid=113068423500393235248&rtpof=true&sd=true
- MBHCatalogs github for code:
  https://github.com/mbonetti90/MBHCatalogs
- MBHCatalogs google drive for simulation data:
  https://drive.google.com/drive/folders/1rW-cdOrMglfp2w72X0jmdu1G5COb7RDk?usp=sharing
- MBHCatalogs general folder
  https://drive.google.com/drive/folders/1YQYw-Km-N5b5jjeHYQ4Ls6LTbXxKh9VC?usp=sharing
Authors
-------
- Matteo Bonetti : matteo.bonetti@unimib.it
- Luke Zoltan Kelley : lzkelley@berkeley.edu

- Structure
-------
This script has four main parts:
    - part A [MODIFY]       : functions that need input from users to collect information about the specific models.
                            Please carefully read the docstring of the function 'input_data()' that you find below.
    - part B [DO NOT MODIFY]: data structure for the hdf5 file
    - part C [DO NOT MODIFY]: functions to produce and validate hdf5 files. 
    - part D [DO NOT MODIFY]: main routine.

----
- Validation
    - Make sure numbers of elements all match
    - Even when not doing units checks, make sure values are sane (e.g. all positive, nonzero; etc) - DONE

---
- Simple usage:

python3 hdf5_catalog_v2.1.py -W -F the_name_of_the_file_to_be_produced

type -h for immediate help


"""
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

# relevant packages, DO NOT REMOVE THEM
import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import h5py
import numpy as np
import argparse
import sys

__VERSION__ = '2.1' # catalog version
DEBUG = False
DEF_FILENAME = "/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/Astrid"
M_DIR = "/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/"
ALL_DIR = "/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/"
np.random.seed(5)

############################################################################################
############################################################################################
############################################################################################
#### PART A: CHANGES NEEDED FROM USERS #####################################################
############################################################################################
############################################################################################
############################################################################################

def input_data():
    '''
    This function collects the necessary information to produce the MBH catalog.

    Users should perform the following actions:
    1) edit the 'metadata' dictionary in this function to provide with the specific information concerning a certain model;
    2) edit the 'MBHB_no_delay' function to provide with the no-delay binary catalog;
    3) [optional] edit the 'MBHB_delay' function to provide with the delay binary catalog. Leave an empty dict if the delay model is not present;
    4) edit the MBH_population function to provide with the complete MBH population at selected redshifts.

    IMPORTANT: If a field does not apply to your specific model (e.g. spatial resolution for EPS SAMs), set it to 'np.nan'.

    NOTE: we require the number density of events. This is a meaninful quantity for EPS SAMs, for cosmological simulations
    and SAMs based on dark matter merger trees this is just 1/V_box 

    All parts that need to be modified by the user are enclosed into starred blocks like this:
    #************************************************#
    #************************************************#
    ...
    ...
    ...
    ...
    #************************************************#
    #************************************************#
        
    All the other functions should be left unchanged.

    Returns: 4 dictionaries -> metadata, allbhs, binaries_nodelay, binaries_delay
    '''

    ###########################################################
    ####### CHANGE THE CODE BELOW #############################
    ###########################################################
    #*********************************************************#
    #*********************************************************#

    metadata = {
        # ---- Header data - identifying information for the dataset
        # REQUIRED
        'SimulationName': 'Astrid',
        'SimulationVersion': '1.0',
        'ModelType': 'Hydro',
        'Version': __VERSION__,
        'Date': str(datetime.datetime.now()),
        'Contributor': ["Nianyi Chen", "Yueying Ni"],
        'Email': ["nianyi.chen7@gmail.com", "yueying.ni@cfa.harvard.edu"],
        'Principal': ["Tiziana Di Matteo"],
        'Reference': ["10.1093/mnras/stac1432", "10.1093/mnras/stac351", "10.1093/mnras/stad3084"],
        # OPTIONAL:
        'Website': ["https://astrid.psc.edu"],

        # ---- Model parameters - metadata specification for simulation(s) used to construct catalog
        # REQUIRED
        'HubbleConstant': 67.74, #km s^-1 Mpc^-1
        'OmegaMatter': 0.3089,
        'OmegaLambda': 0.6911,
        'BoxSize': 369, # cMpc
        'MinBHSeedMass': 4.4e4, #M_sun
        'MinRedshift': 0.2,
        'MaxRedshift': 99,
        'StellarMassResolution': 4.7e5, # M_sun
        'DarkMatterMassResolution': 9.63e6, # M_sun
        'SpatialResolution': 1.5, # kpc
        # OPTIONAL:
        'MinimumDarkMatterHaloMass': 3e8, # M_sun
        'GasMassResolution': 1.87e6, # M_sun
    }
    #*********************************************************#
    #*********************************************************#

    ###########################################################
    #### DO NOT CHANGE FUCNTION CALLS AND RETURN ##############
    ###########################################################

    # REQUIRED
    # ---- Binary mergers at time of galaxy merger (no-delay) ----

    binaries_nodelay = MBHB_no_delay(metadata)


    # OPTIONAL (can be empty dictionary):
    # ---- Binary mergers including a model of delay times ----

    binaries_delay = MBHB_delay(metadata)


    # REQUIRED:
    # ---- All black holes and host galaxies at target redshifts

    allbhs = MBH_population(metadata)

    #########################

    return metadata, allbhs, binaries_nodelay, binaries_delay

############################################################################################

def MBHB_no_delay(metadata):

    '''
    Function to collect properties of binaries assuming no-delay models.

    THE CODE BELOW PRODUCES FAKE BINARIES PROPERTIES, 
    PLEASE REPLACE IT WITH CUSTOMARY CODE TO COLLECT BINARIES FROM YOUR MODEL.

    # Required fields
        N_binaries: number of binaries 
        m1: primary mass [M_sun]
        m2: secondary mass [M_sun]
        z: redshift of binary merger [None]
        sepa: separation at merger [proper kpc]
        W: number density [Mpc^-3]
        mstar: host galaxy stellar mass at merger or just after it [M_sun]
        zgal: redshift of galaxy stellar mass [None]
        R50: half-mass radius or effective radius [proper kpc]
    
    # optional fields
        mdm: dark matter mass at merger or just after it [M_sun]

    Returns: dictionary 
    '''

    ###########################################################
    ####### CHANGE THE CODE BELOW #############################
    ###########################################################
    #*********************************************************#
    #*********************************************************#

    # generate fake binary properties
    mdata = np.load(M_DIR + "lisa_mcat_v2.1_z0.57.npy")

    N_binaries = len(mdata)
    m1 = mdata["m1"]
    m2 = mdata["m2"]
    m1, m2 = np.max([m1, m2], axis=0), np.min([m1, m2], axis=0)
    
    z = mdata["z"]
    sepa = mdata["dr"]  # in kpc unit

    # generate fake binary-host galaxy-remnant properties
    mstar = mdata["mstot"]
    mdm = -np.ones(len(mdata["mstot"]))
    zgal = mdata["zsnap"]
    print("zgal:", mdata["zsnap"][:10])
    R50 = np.random.normal(loc=3, scale=0.1, size=N_binaries)

    # number density
    W = 1/metadata["BoxSize"]**3*np.ones(len(m1))

    # FILL METADATA INFO
    # total number of merged binaries assuming no delays
    metadata['NumberBinaries'] = N_binaries

    # provide an explanation of the merger criterion, modify the string
    metadata['NoDelay'] = (
        "mergers occur when two MBH particles come within 2*gravitational softening length of "
        "eachother, and the kinetic energy of the pair is less than the gravitational "
        "potential energy between them.")
    
    metadata['CommentsNoDelay'] = (
        "The simulation contains subgrid dynamical friction model, so separation should be used as the initial separation for post-processing instead of galaxy sizes")

    #**********************************************************#
    #**********************************************************#
    ############################################################
    # DO NOT CHANGE DICTIONARY AND RETURN ######################
    ############################################################

    # collect data in a dictionary
    binaries_nodelay = {
        "BlackHoles": {
            'PrimaryMass': m1,
            'SecondaryMass': m2,
            'Redshift': z,
            'Separation': sepa,
            'NumberDensity': W,
        },
        "Galaxies": {
            'RemnantStellarMass': mstar,
            'RemnantRedshift': zgal,
            'RemnantR50': R50,
            # OPTIONAL:
            'RemnantDarkMatterMass': mdm,
        }
    }

    return binaries_nodelay

#########################
def MBHB_delay(metadata):

    '''
    Function to collect properties of binaries assuming delay models.

    THE CODE BELOW PRODUCES FAKE BINARIES PROPERTIES, 
    PLEASE REPLACE IT WITH CUSTOMARY CODE TO COLLECT BINARIES FROM YOUR MODEL.
    IF YOU DO NOT PROVIDE ANY DELAY MODEL JUST SET 'delay_model' TO 'False' AND DELETE THE CODE, AN EMPTY DICT WILL BE SET.

    Returns:  
    # Required fields
        N_binaries_delay: number of binaries 
        m1_delay: primary mass [M_sun]
        m2_delay: secondary mass [M_sun]
        z_delay: redshift of binary merger [None]
        sepa_delay: separation at merger [proper kpc]
        W: number density [Mpc^-3]
        mstar_delay: host galaxy stellar mass at merger or just after it [M_sun]
        zgal_delay: redshift of galaxy stellar mass [None]
        R50: half-mass radius or effective radius [proper kpc]
    
    # optional fields
        mdm_delay: dark matter mass at merger or just after it [M_sun]
    '''

    ###########################################################
    ####### CHANGE THE CODE BELOW #############################
    ###########################################################
    #*********************************************************#
    #*********************************************************#

    # set to true if data for a model with delays are available
    delay_model = False

    # generate fake binary properties
    N_binaries = 1000
    m1 = np.random.normal(loc=1e7, scale=1e6, size=N_binaries)
    m2 = np.random.normal(loc=1e7, scale=1e6, size=N_binaries)
    m1, m2 = np.max([m1, m2], axis=0), np.min([m1, m2], axis=0)
    z = 0.01 + np.random.uniform(0, 10, size=N_binaries)
    sepa = 10.0 ** np.random.uniform(2.0, 4.0, size=N_binaries)

    # generate fake binary-host galaxy-remnant properties
    mstar = np.random.normal(loc=1e10, scale=1e8, size=N_binaries)
    mdm = mstar * np.maximum(1.0, np.random.normal(loc=10.0, scale=1.0, size=N_binaries))
    zgal = z - np.random.uniform(0, 0.001, size=N_binaries)
    R50 = np.random.normal(loc=3, scale=1, size=N_binaries)

    # choose some amount of delay in redshift
    delta_redz = np.random.uniform(0.01, 1.0, N_binaries)
    z_delay = z - delta_redz
    zgal_delay = zgal - delta_redz
    # valid binaries must coalesce before redshift zero
    sel_valid = (z_delay > 0.0) & (zgal_delay > 0.0)
    m1_delay = m1[sel_valid]
    m2_delay = m2[sel_valid]
    z_delay = z_delay[sel_valid]
    zgal_delay = zgal_delay[sel_valid]
    sepa_delay = sepa[sel_valid] * (10.0 ** np.random.normal(loc=0.02, scale=0.05))
    # make sure delayed separation is smaller than no-delay separation
    sepa_delay = np.minimum(sepa_delay, sepa[sel_valid]*0.98)

    # find number of coalesced binaries in the delay model
    N_binaries_delay = int(m1_delay.size)

    # define new arrays with delayed binaries
    m1_delay = m1_delay * (1.0 + np.random.normal(loc=0.1, scale=0.1, size=N_binaries_delay))
    m2_delay = m2_delay * (1.0 + np.random.normal(loc=0.2, scale=0.1, size=N_binaries_delay))
    mstar_delay = mstar[sel_valid] * (1.0 + np.random.normal(loc=0.15, scale=0.1, size=N_binaries_delay))
    mdm_delay = mstar[sel_valid] * (1.0 + np.random.normal(loc=0.15, scale=0.1, size=N_binaries_delay))
    # primary and secondary may have switched
    m1_delay, m2_delay = np.max([m1_delay, m2_delay], axis=0), np.min([m1_delay, m2_delay], axis=0)

    # number density
    W = 1/metadata["BoxSize"]**3*np.ones(len(m1_delay))

    # FILL METADATA INFO
    # total number of merged binaries assuming delays
    metadata["NumberBinariesDelay"] = N_binaries_delay

    # provide an explanation of the merger criterion, modify the string
    metadata["Delay"] = (
        "Binaries are evolved in post-processing using a semi-analytic model for dynamical "
        "friction.  Binaries are considered to have merged once they become gravitationally bound."
    )

    metadata['CommentsDelay'] = (
        "Any additional information you consider relevant for any clarification, i.e."
        "model special features, recipe to deal with MBH evolution etc.")

    #**********************************************************#
    #**********************************************************#
    ############################################################
    # DO NOT CHANGE DICTIONARY AND RETURN ######################
    ############################################################

    if delay_model:
        # collect data in a dictionary
        binaries_delay = {
            "BlackHoles": {
                'PrimaryMass': m1_delay,
                'SecondaryMass': m2_delay,
                'Redshift': z_delay,
                'Separation': sepa_delay,
                'NumberDensity': W,
            },
            "Galaxies": {
                'RemnantStellarMass': mstar_delay,
                'RemnantRedshift': zgal_delay,
                'RemnantR50': -np.ones(len(R50)),
                # OPTIONAL:
                'RemnantDarkMatterMass': mdm_delay,
            }
        }

    # overwrite with an emtpy dict if no model with delays is present
    else:
        binaries_delay = {}

    return binaries_delay

#########################
def MBH_population(metadata):

    '''
    Function to collect properties of all BHs at selected redshift [TARGET_REDSHIFTS array].
    See also https://docs.google.com/spreadsheets/d/1Zz-A5QcKHdPabpSmkKcs9ERzeAA0sqM_/edit#gid=1167697008

    THE CODE BELOW PRODUCES FAKE BHs PROPERTIES, 
    PLEASE REPLACE IT WITH CUSTOMARY CODE TO COLLECT BHs FROM YOUR MODEL.

    # Required fields at selected target redshifts:
        mass: BH mass [M_sun]
        mdot: accretion rate [M_sun/yr]
        radeff: radiative efficiency [-]
        W: number density [Mpc^-3]
        mstar: host galaxy stellar mass [M_sun]
        mdm: dark matter mass at merger [M_sun]


    Returns: dictionary 
    '''

    # data structure to contain collected data
    allbhs = {}
    allbh_redshifts = []
    num_allbhs = 0

    ###########################################################
    ####### CHANGE THE CODE BELOW #############################
    ###########################################################
    #*********************************************************#
    #*********************************************************#
    
    min_redshift = np.inf
    for tarz in TARGET_REDSHIFTS:
        if tarz == 5:
            data = np.load(ALL_DIR + "lisa_mbhcat_single_v2.1_z5.npy")
            print("loading from:", ALL_DIR + "lisa_mbhcat_single_v2.1_z5.npy")
        elif tarz == 2:
            data = np.load(ALL_DIR + "lisa_mbhcat_single_v2.1_z2.npy")
            print("loading from:", ALL_DIR + "lisa_mbhcat_single_v2.1_z2.npy")
        elif tarz == 1:
            data = np.load(ALL_DIR + "lisa_mbhcat_single_v2.1_z1.npy")
            print("loading from:", ALL_DIR + "lisa_mbhcat_single_v2.1_z1.npy")
        else:
            continue
        myz = tarz
        # store actual simulation redshifts
        allbh_redshifts.append(myz)
        myz = np.maximum(myz, 1.0e-4)
        # generate a key for this redshift (i.e. string specification in consistent format)
        zkey = f"Z{myz:011.8f}"
        # number of BHs in this snapshot
        numbhs = len(data)
        # store the number of BHs in the lowest-redshift snapshot being stored
        if myz < min_redshift:
            min_redshift = myz
            num_allbhs = numbhs

        # generate fake BH data
        mass = data["bhmass"]
        mdot = data["mdot"]
        radEff = 0.1*np.ones(len(mdot)) # example with fixed rad efficiency for all BHs
        # generate fake galaxy data
        mstar = data["mstot"]
        mdm = data["mdm"]

        # number density
        W = 1/metadata["BoxSize"]**3*np.ones(len(mass))

        # put data for this redshift/snapshot into nested dictionaries
        temp_data = {
            "BlackHoles": {
                "Mass": mass,
                "Mdot": mdot,
                "RadEff": radEff,
                "NumberDensity": W,
            },
            "Galaxies": {
                "StellarMass": mstar,
                "DarkMatterMass": mdm,
            },
        }
        # add to combined allbhs data
        allbhs[zkey] = temp_data
        print(zkey)

    metadata['CommentsAllBHs'] = (
        "Any additional information you consider relevant for any clarification, i.e."
        "model special features, recipe to deal with MBH evolution etc.")

    #*********************************************************#
    #*********************************************************#
    ###########################################################
    # DO NOT CHANGE METADATA ASSIGNMENT AND RETURN ############
    ###########################################################

    metadata['NumberAllBHs'] = num_allbhs
    metadata['RedshiftsAllBHs'] = allbh_redshifts

    return allbhs

############################################################################################
############################################################################################
############################################################################################
#### PART B: DATA STRUCTURE, DO NOT MODIFY #################################################
############################################################################################
############################################################################################
############################################################################################

'''
The following dictionaries contain as keys the names of the fields that populate the hdf5 files. 
The values associated with keys, either a list or a scalar, are used to validate the catalog and do not
enter the catalog itself (except for units).
When values are lists, those are formed by three elements:
- string or None: the units of the quantity (check with https://docs.google.com/spreadsheets/d/1Zz-A5QcKHdPabpSmkKcs9ERzeAA0sqM_/edit#gid=1167697008)
- float: lower bound (can be None meaning no bound)
- float: upper bound (can be None meaning no bound)
'''

# binaries
DATA_BINARIES = {
    "BlackHoles": {
        "PrimaryMass": ["M_sun", 1.0e2, 1.0e12],
        "SecondaryMass": ["M_sun", 1.0e2, 1.0e12],
        "Redshift": [None, 0.0, 100.0],
        "Separation": ["kpc", 0.0, 1.0e4],
        "NumberDensity": ["Mpc^-3", 1e-30, 1e10],
    },
    "Galaxies": {
        "RemnantStellarMass": ["M_sun", 1.0e4, 1.0e14],
        "RemnantDarkMatterMass": ["M_sun", 1.0e4, 1.0e16],
        "RemnantRedshift": [None, 0.0, 40.0],
        "RemnantR50": [None, 1e-3, 1e3],
    },
}

DATA_BINARIES_OPTIONAL = ["RemnantDarkMatterMass","RemnantR50"]


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

DATA_ALLBHS_OPTIONAL = ["DarkMatterMass"]

TARGET_REDSHIFTS = [0, 1, 2, 5, 10, 15, 20]


# this dictionary with fields for the hdf5 header
METADATA = {
    # Header
    'SimulationName': None,
    'SimulationVersion': None,
    'NumberBinaries': None,
    'NumberBinariesDelay': None,
    'NumberAllBHs': None,
    'RedshiftsAllBHs': [None, 0, 100.0],
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
    'OmegaMatter': [None,0,1],
    'OmegaLambda': [None,0,1],
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

METADATA_OPTIONAL = ["Website", "GasMassResolution",'MinimumDarkMatterHaloMass',"NumberBinariesDelay","Delay","CommentsNoDelay","CommentsDelay","CommentsAllBHs"]

ERROR_LIST = []
WARN_LIST = []

############################################################################################
############################################################################################
############################################################################################
#### PART C: NO CHANGES NEEDED FROM USERS, DO NOT MODIFY CODE BELOW ########################
############################################################################################
############################################################################################
############################################################################################

def write_hdf5(fname, metadata, allbhs, binaries_nodelay, binaries_delay):
    '''
    This function creates a hdf5 file containing the properties of the catalog. The file has the following structure:
    - Header group

    - Binaries group
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
                - RemnantR50

        - Delay group
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


        # ---- write binaries data to "Binaries/" ----
        binaries_group = h5out.create_group("Binaries")

        nodelay_group = binaries_group.create_group("NoDelay")
        # iterate over "BlackHoles" and "Galaxies"
        for group_key, group_vals in binaries_nodelay.items():
            # e.g. "Binaries/NoDelay/BlackHoles"
            group = nodelay_group.create_group(group_key)
            # iterate over each element of these groups
            # e.g. "PrimaryMass" ...  and "RemanntStellarMass" ...
            for key, vals in group_vals.items():
                print("key0", key)
                if DEBUG:
                    print(f"writing {group.name}[{key}] = {type(vals)}, {vals.shape}, {vals.dtype}")
                # e.g. "Binaries/NoDelay/BlackHoles/PrimaryMass"
                tmp = group.create_dataset(key, data=vals)
                print("key", key)
                tmp.attrs.create('units',data=str(DATA_BINARIES[group_key][key][0]))
                print("key2", key)
            print("flag1")
        
        if (binaries_delay is not None) and (len(binaries_delay.keys()) > 0):
            delay_group = binaries_group.create_group("Delay")
            # iterate over "BlackHoles" and "Galaxies"
            for group_key, group_vals in binaries_delay.items():
                # e.g. "Binaries/Delay/BlackHoles"
                group = delay_group.create_group(group_key)
                # iterate over each element of these groups
                # e.g. "PrimaryMass" ...  and "RemanntStellarMass" ...
                for key, vals in group_vals.items():
                    if DEBUG:
                        print(f"writing {group.name}[{key}] = {type(vals)}, {vals.shape}, {vals.dtype}")
                    # e.g. "Binaries/Delay/BlackHoles/PrimaryMass"
                    tmp = group.create_dataset(key, data=vals)
                    tmp.attrs.create('units',data=str(DATA_BINARIES[group_key][key][0]))

        # ---- write All BHs data to "AllBHs/" ----
        allbhs_group = h5out.create_group("AllBHs")
        for zkey, zdata in allbhs.items():
            # e.g. "AllBHs/Z00.00215210"
            redz_group = allbhs_group.create_group(zkey)
            # iterate over "BlackHoles" and "Galaxies"
            for group_key, group_data in zdata.items():
                # this is ""AllBHs/Z00.00215210/BlackHoles"  or  "AllBHs/Z00.00215210/Galaxies"
                group = redz_group.create_group(group_key)
                # iterate over each element of these groups
                # e.g. "Mass" ...  and "RemanntStellarMass" ...
                for key, vals in group_data.items():
                    if DEBUG:
                        print(f"writing {group.name}[{key}] = {type(vals)}, {vals.shape}, {vals.dtype}")
                    # e.g. "AllBHs/Z00.00215210/BlackHoles/Mass"
                    tmp = group.create_dataset(key, data=vals)
                    tmp.attrs.create('units',data=str(DATA_ALLBHS[group_key][key][0]))

    if DEBUG:
        print("\nDone")
        print("###########################################\n")

    return

##############################################################################

def read_hdf5(name, h5in, indent=0):
    """Recursively print out the contents of the given HDF5 element."""

    # print with indentation to show heirarchy of elements
    def my_print(vals, indent):
        pref = " " * (indent+1) * 4
        print(f"{pref}{vals}")

    # if this is a group or other container, print its name
    try:
        # this will raise an error if this is NOT a group/container
        h5in.items()
        my_print(f"{name}", indent-1)
    except AttributeError:
        pass

    # if there is header (`attrs`) information, print it
    header_keys = h5in.attrs.keys()
    if len(header_keys) > 0:
        my_print("attrs", indent)
        for key, vals in h5in.attrs.items():
            my_print(f"{key}: {vals}", indent+1)

    # try to recursively descend into groups below this level
    try:
        for key, vals in h5in.items():
            read_hdf5(key, vals, indent+1)

    # if we can't descend further, this is a dataset, print its basic information
    except AttributeError:
        my_print(f"{name:20s} : {np.shape(h5in)} {h5in.dtype} [{np.min(h5in):+5.2e}, {np.max(h5in):+5.2e}]", indent)
        print("\n")

    return

##############################################################################

def plot_hdf5(name, h5in, choice):

    if choice == "AllBHs":

        data = h5in[choice]

        for k1 in data.keys():
            for k2 in data[k1].keys():
                for k3 in data[k1][k2].keys():

                    if k3 in ["Mass","StellarMass","DarkMatterMass"]:
                        dd = np.log10(data[k1][k2][k3])
                        ll = r"$\log_{10}$ "+k3
                    else:
                        dd = data[k1][k2][k3]
                        ll = k3

                    fig, ax=plt.subplots()
                    ax.set_title(f"{choice}: {k1}")
                    ax.hist(dd,bins=50,range=(np.min(dd),np.max(dd)),histtype='step')
                    ax.set_xlabel(ll)

    if choice == "NoDelay":

        data = h5in["Binaries"][choice]

        for k1 in data.keys():
            for k2 in data[k1].keys():

                    if k2 in ["PrimaryMass","SecondaryMass","RemnantStellarMass","RemnantDarkMatterMass"]:
                        dd = np.log10(data[k1][k2])
                        ll = r"$\log_{10}$ "+k2
                    else:
                        dd = data[k1][k2]
                        ll = k2

                    fig, ax=plt.subplots()
                    ax.set_title(f"{choice}: {k1}")
                    ax.hist(dd,bins=50,range=(np.min(dd),np.max(dd)),histtype='step')
                    ax.set_xlabel(ll)

    if choice == "Delay":
        data = h5in["Binaries"][choice]

        for k1 in data.keys():
            for k2 in data[k1].keys():

                    if k2 in ["PrimaryMass","SecondaryMass","RemnantStellarMass","RemnantDarkMatterMass"]:
                        dd = np.log10(data[k1][k2])
                        ll = r"$\log_{10}$ "+k2
                    else:
                        dd = data[k1][k2]
                        ll = k2

                    fig, ax=plt.subplots()
                    ax.set_title(f"{choice}: {k1}")
                    ax.hist(dd,bins=50,range=(np.min(dd),np.max(dd)),histtype='step')
                    ax.set_xlabel(ll)
    return

##############################################################################

def get_filename(args):

    fname = args.filename.strip().strip(".hdf5")
    fname = fname + f"_v{__VERSION__}.hdf5"
    fname = Path(fname).resolve()
    if not fname.parent.is_dir():
        raise FileNotFoundError(f"Path to file '{fname.parent}' does not exist!")

    return fname

##############################################################################

def validate(h5in):

    import re

    error_flag = 0
    warn_flag = 0

    def error(msg):
        nonlocal error_flag
        error_flag += 1
        print(f"ERROR  : {msg} !!")
        ERROR_LIST.append(msg)

    def warn(msg):
        nonlocal warn_flag
        warn_flag += 1
        print(f"WARNING: {msg}")
        WARN_LIST.append(msg)

    def note(msg):
        print(f"{msg}")

    def check_units_bounds(key, vals, check):
        """If `check` is a (3,) iterable, make sure all values are within the bounds specified.

        `check`:  [units string, lower-bound float, upper-bound float]
            e.g. ["M_sun", 1.0e2, 1.0e10]  means solar-mass units, between 1e2 and 1e10.
            The lower-bound and upper-bound can also be given as `None` to signify no bound.

        """

        if np.size(check) != 3:
            return
                
        unit, lo, hi = check
        good_lo = (lo is None) or np.all(np.greater_equal(vals, lo))
        good_hi = (hi is None) or np.all(np.less_equal(vals, hi))
        if not good_lo or not good_hi:
            nan_val = False
            aa = np.isnan(vals)
            if len(aa[aa==True]) > 0:
                nan_val = True
            print("\n")
            warn(f"Units and bounds: {key} values are outside of expectations [{lo:.2e}, {hi:.2e}] units: {unit}. Check low limit passed:{good_lo}, check high limit passed:{good_hi}. Check for nan:{nan_val}. If nan are present make sure this was intentional!")

        return

    def check_expected_exist(namespace, specs, name, optional=[]):
        """`specs` specifies the key-value pairs of `namespace`.  Make sure required elements are defined.
        """
        # make sure required elements are present
        for key, check_vals in specs.items():
            required = (key not in optional)
            exists = (key in namespace.keys())
            if DEBUG:
                print(f"check: {key} in {name} :: required={required}, exists={exists}")
            if exists:
                if DEBUG:
                    key2 = str(namespace)+key 
                else:
                    key2 = key
                check_units_bounds(key2, namespace[key], check_vals)
            elif required:
                error(f"{name} key `{key}` not in {namespace}")

        # see if unexpected keys are present
        for key, vals in namespace.items():
            if (key not in specs.keys()):
                note(f"found unexpected `{namespace}`: {key}={vals}")

    # ---- validate meta-data (stored in Header as `attrs`)

    # make sure `METADATA` specifications are defined in hdf5 Header
    check_expected_exist(h5in['Header'].attrs, METADATA, "metadata", METADATA_OPTIONAL)

    # ---- validate `Binaries` data

    if "Binaries" in h5in:
        binaries = h5in['Binaries']
        delay_types = ["NoDelay", "Delay"]
        is_required = [True, False]
        # iterate over "Binaries/NoDelay" and "Binaries/Delay"
        for _delay, reqd in zip(delay_types, is_required):
            if _delay in binaries:
                dgroup = binaries[_delay]

                # iterate over "[No]Delay/BlackHoles" and "[No]Delay/Galaxies" along with their specifications
                for bgroup_key, bgroup_vals in DATA_BINARIES.items():

                    check_expected_exist(
                        dgroup[bgroup_key],  # make sure this group (e.g. "[No]Delay/BlackHoles")
                        bgroup_vals,         # contains these values
                        f"Binaries/{_delay}/{bgroup_key}",
                        DATA_BINARIES_OPTIONAL
                    )

            else:
                msg = f"Key `{_delay}` not found in `Binaries`"
                if reqd:
                    error(msg)
                else:
                    note(msg)

        # check for negative or nan values
        for _delay in delay_types:
            dgroup = binaries[_delay]
            for _key in dgroup.keys():
                ddgroup = dgroup[_key]
                for key in ddgroup.keys():
                    tmp = ddgroup[key][:]
                    if len(tmp[tmp<0]) > 0 or len(tmp[tmp==np.nan]) > 0:
                        msg = f"{_delay}/{_key}/{key} -> Negative or nan value(s)"
                        error(msg)

        # check that primary mass is larger than secondary mass
        for _delay in delay_types:
            dgroup = binaries[_delay]
            m1 = dgroup["BlackHoles/PrimaryMass"][:]
            m2 = dgroup["BlackHoles/SecondaryMass"][:]
            z = dgroup["BlackHoles/Redshift"][:]

            ch = np.greater_equal(m1,m2)
            if len(ch[ch==False]) > 0:
                msg = f"Binaries/{_delay}/BlackHoles/PrimaryMass must be >= than {_delay}/BlackHoles/SecondaryMass"
                error(msg)

            if len(z[z<0]) > 0:
                msg = f"Binaries/{_delay}/BlackHoles/Redshift must be positive"
                error(msg)

    else:
        error("`Binaries` entry not found in file!")

    # ---- validate AllBHs data

    if "AllBHs" in h5in:
        allbhs = h5in['AllBHs']

        zkeys = allbhs.keys()
        target_found = np.zeros(len(TARGET_REDSHIFTS), dtype=bool)
        for ii, zk in enumerate(zkeys):
            if DEBUG:
                print(f"`AllBHs` entry {ii}: '{zk}' (length={len(zk)})")

            # make sure that redshift entry looks correctly formatted, with the correct precision
            try:
                assert (len(zk) >= 12)
                match = re.match(r'Z(\d{2}\.\d{8})', zk)
                redz = float(match.group(1))
            except Exception:
                error(f"Entry in `AllBHs` (`{zk}`) is unexpected.  Expected 'Z<%011.8f>' e.g. 'Z01.99219512'")
                continue

            # check that redshift entry is near one of the requested target redshifts
            found = np.isclose(redz, TARGET_REDSHIFTS, rtol=0.1, atol=0.5)
            if not np.any(found):
                warn(f"`AllBHs` entry '{zk}' ({redz}) is not close to target redshifts {TARGET_REDSHIFTS}")
            # mark that this target redshift was found
            target_found[found] = True

            zgroup = allbhs[zk]
            # iterate over "BlackHoles" and "Galaxies" and check each of their specifications
            for allbhs_key, allbhs_vals in DATA_ALLBHS.items():
                check_expected_exist(
                    zgroup[allbhs_key],  # make sure this group (e.g. "AllBHs/Z01.99219512/BlackHoles")
                    allbhs_vals,         # contains these values
                    f"AllBHs/{zk}/{allbhs_key}",
                    DATA_ALLBHS_OPTIONAL,
                )
                
            for allbhs_key, allbhs_vals in zgroup.items():
                for key, vals in zgroup[allbhs_key].items():
                    tmp = vals[:]
                    if len(tmp[tmp<0]) > 0 or len(tmp[tmp==np.nan]) > 0:
                        msg = f"AllBHs/{zk}/{allbhs_key}/{key} -> Negative or nan value(s)"
                        error(msg)

        # make sure all of the target redshifts have been found
        if not np.all(target_found):
            not_found = np.array(TARGET_REDSHIFTS)[~target_found]
            warn(f"Missing {len(not_found)} target redshifts: {not_found}")

    else:
        error("`AllBHs` entry not found in file!")


    if error_flag:
        print("\n\n######################")
        print("ERROR LIST:")
        for er in ERROR_LIST:
            print(er)
        print("######################\n\n")
        raise RuntimeError(f"File {h5in} is invalid!  Found {error_flag} critical errors.")
    
    if warn_flag:
        print("\n\n######################")
        print("WARNING LIST:")
        for wa in WARN_LIST:
            print(wa)
        print("######################\n\n")

    return


############################################################################################
############################################################################################
############################################################################################
#### PART D: MAIN ROUTINE, NO CHANGES NEEDED FROM USERS, DO NOT MODIFY CODE BELOW ##########
############################################################################################
############################################################################################
############################################################################################

if __name__ == "__main__":
    # Parser instructions
    parser = argparse.ArgumentParser(description='Catalog file creator')
    parser.add_argument('-D', '--debug', action="store_true", default=False, help='Run in debug mode. This enables many prints during file writing (default False)')
    parser.add_argument('-W', '--write', action="store_true", default=False, help='Write the hdf5 file (default False)')
    parser.add_argument('-R', '--read', action="store_true", default=False, help='Read the hdf5 file (default False)')
    parser.add_argument('-P', '--plot', default=False,  choices=["NoDelay","Delay","AllBHs"], help='Plot histograms of datasets in the hdf5 file. Choose among NoDelay, Delay or AllBHs.')
    parser.add_argument('-F', '--filename', default=DEF_FILENAME, help='Name of the hdf5 file (default example-sim)')
    parser.add_argument('-V', '--validate', action="store_true", default=False, help='Check HDF5 file for consistency and completeness.')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    if args.debug:
        DEBUG = True

    # name of the hdf5 file to be produced
    fname = get_filename(args)

    # write the file
    if args.write:
        print("\nLoading input data")
        metadata, allbhs, binaries_nodelay, binaries_delay = input_data()
        print(f"\nWriting file {fname}")
        write_hdf5(fname, metadata, allbhs, binaries_nodelay, binaries_delay)
        print("\nFile has been written")

    if not fname.is_file():
        raise FileNotFoundError(f"File {fname} does not exist!")

    # read the file to check data
    if args.read:
        print(f"\nReading file",end=" ")
        with h5py.File(fname, 'r') as h5in:
            read_hdf5(fname, h5in)
        print("\nFile successfully read")

    # verify file contents
    if args.validate:
        print(f"\nValidating file {fname}")
        with h5py.File(fname, 'r') as h5in:
            validate(h5in)
        print("\nValidation completed")

    # make some plots
    if args.plot:
        with h5py.File(fname, 'r') as h5in:
            plot_hdf5(fname, h5in, args.plot)
        plt.show()
        print("\nPlotting completed")


    print("")
