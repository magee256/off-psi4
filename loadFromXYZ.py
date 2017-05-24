#!/usr/bin/env python


## Write new coordinates from an XYZ file to an existing mol2
## file with old coordinates. Helpful for VMD optimization and
## keeping consistent atom numbering for RMSD comparisons.

## Usage: python loadFromXYZ.py -m initial.mol2 -x final.xyz -o final.mol2

import os, sys
import openeye.oechem as oechem
import argparse


### ------------------- Script -------------------

def main(**kwargs):
    outfn = opt['fout']
    
    # Open input files.
    mifs = oechem.oemolistream()
    if not mifs.open(opt['fmol2']):
        oechem.OEThrow.Warning("Unable to open %s for reading" % opt['fmol2'])
        return
    xifs = oechem.oemolistream()
    if not xifs.open(opt['fxyz']):
        oechem.OEThrow.Warning("Unable to open %s for reading" % opt['fxyz'])
        return
    mmol = mifs.GetOEMols().next()
    xmol = xifs.GetOEMols().next()

    mmol.SetCoords(xmol.GetCoords())
    
    ofs = oechem.oemolostream()
    if not ofs.open(outfn):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfn)
    oechem.OEWriteConstMolecule(ofs, mmol)
    ofs.close()
    mifs.close()
    xifs.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--fmol2",
        help="Reference mol2 file which new coordinates will be loaded into.")
    parser.add_argument("-x", "--fxyz",
        help="XYZ file with new coordinates to load into mol2 file.")
    parser.add_argument("-o", "--fout",
            help="Name of the output mol2 file.")

    args = parser.parse_args()
    opt = vars(args)

    main(**opt)

