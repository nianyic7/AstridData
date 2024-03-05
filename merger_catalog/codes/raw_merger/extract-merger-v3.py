"""
extract merger event from bhdetail-reduce files
Usage: python extract-merger-v3.py $snap

An parallelized version of extract-merger-v2.py,
contains additional information about mergers:
positions and accretion rates just before the merger.
"""
import numpy as np
from bigfile import BigFile
import sys
import glob


hh = 0.6774
mdot_msun_yr = 1e10/978/1e6

class ChunkData:
    def __init__(self, bf, ichunk):
        self.bf = bf
        self.ichunk = ichunk
        self.acBHMass_ds = None
        self.swallowed_ds = None
        self.swallowID_ds = None
        self.BHID_ds = None
        self.redshift_ds = None
        self.BHpos_ds = None
        self.BHvel_ds = None
        self.BHMdot_ds = None
        self.BHMass_ds = None
        self.zstart = None
        self.zend = None

    def load(self):
        self.acBHMass_ds = bf['acBHMass']
        self.swallowed_ds = bf['swallowed']
        self.swallowID_ds = bf['swallowID']
        self.BHID_ds = bf['ID']
        self.redshift_ds = bf['redshift']
        self.BHpos_ds = bf['Coordinates']
        self.BHvel_ds = bf['Velocity']
        self.BHMdot_ds = bf['AccretionRate']
        self.BHMass_ds = bf['BH_Mass']
        self.zstart = bf.attrs['zstart']
        self.zend = bf.attrs['zend']
