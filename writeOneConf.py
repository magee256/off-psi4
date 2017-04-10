#!/usr/bin/env python


## This script writes out a single conformer of some mol from
## an SDF file having multiple molecules, and multiple conformers
## of each molecule. This can be used, for example, to isolate
## specific conformers that lead to a high RMSD energy value 
## for a particular molecule.

## Usage: python writeOneConf.py -f inputfile.sdf -t molName -s sdtag -v tagvalue -x filesuffix

import os, sys
import openeye.oechem as oechem
import argparse


### ------------------- Script -------------------

def main(**kwargs):
    outfn = os.path.splitext(opt['infn'])[0]+'_'+opt['suffix']+'.mol2'
    success = False
    
    ### Read in .sdf file and distinguish each molecule's conformers
    ifs = oechem.oemolistream()
    ifs.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifs.open(opt['infn']):
        oechem.OEThrow.Warning("Unable to open %s for reading" % opt['infn'])
        return
    
    for mol in ifs.GetOEMols():
        if mol.GetTitle() == opt['title']:
            for i, conf in enumerate( mol.GetConfs()):
                if oechem.OEGetSDData(conf, opt['sdtag']) == opt['value']:
                    #print oechem.OEGetSDData(conf, opt['sdtag'])
                    #print opt['value']
                    ofs = oechem.oemolostream()
                    #if os.path.exists(outfn):
                    #    print "Output .sdf file already exists. Exiting.\n"
                    #    return
                    if not ofs.open(outfn):
                        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfn)
                    oechem.OEWriteConstMolecule(ofs, conf)
                    ofs.close()
                    success = True
    if not success:
        print("\n** Found no confs matching your criteria. **")
    ifs.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--infn",
        help="Name of input SDF file from which to write out the conformer.")
    parser.add_argument("-t", "--title",
        help="Molecule name in the file.")
    parser.add_argument("-s", "--sdtag",
        help="SD tag to search for conformer.")
    parser.add_argument("-v", "--value",
        help="Value of the SD tag to write out that conformer.")
    parser.add_argument("-x", "--suffix",
        help="Suffix appened to input fn when writing out this conf.")

    args = parser.parse_args()
    opt = vars(args)

    ### Check that input file exists.
    if not os.path.exists(opt['infn']):
        raise parser.error("Input file %s does not exist. Try again." % opt['infn'])

    main(**opt)

