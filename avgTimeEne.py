#!/usr/bin/env python

# TODO
# - remove GetSDList here and source it from scratch.py

# Notes on error checking:
# - Getting relative energies- ValueError: could not convert string to float
#   Make sure you specify --efromopt if not an SPE file.

import os
import openeye.oechem as oechem
import numpy as np
import argparse

### ------------------- Functions -------------------


def GetSDList(Mol, Property, Method=None, Basisset=None):
    """
    Get list of specified SD tag for all confs in Mol.

    Parameters
    ----------
    Mol:        OEChem molecule with all of its conformers
    Property    string description of property of interest
        options implemented: "QM opt energy" "MM opt energy"
    Method:     optional string, e.g. 'mp2'
    Basisset:   optional string, e.g. '6-31+G(d)'

    Returns
    -------
    sdlist: A 1D N-length list for N conformers with property from SDTag.
    """

    if Property=="QM opt energy":
        taglabel = "QM Psi4 Final Opt. Energy (Har) %s/%s" % (Method, Basisset)

    if Property=="QM spe":
        taglabel = "QM Psi4 Single Pt. Energy (Har) %s/%s" % (Method, Basisset)

    if Property=="MM opt energy":
        taglabel = "MM Szybki Newton Energy"

    if Property=="original index":
        taglabel = "Original omega conformer number"

    if Property=="runtime":
        taglabel = "QM Psi4 Opt. Runtime (sec) %s/%s" % (Method, Basisset)

    if Property=="step":
        taglabel = "QM Psi4 Opt. Steps %s/%s" % (Method, Basisset)

    SDList = []
    for j, conf in enumerate( Mol.GetConfs() ):
        if "DID NOT FINISH" not in oechem.OEGetSDData(conf, "Note on opt.")\
 or Property =="original index":
            SDList.append(oechem.OEGetSDData(conf, taglabel))
        else:
            SDList.append('nan')
    return SDList



def timeAvg(sdfRef, method, basis, steps=False):
    """

    For an SDF file with all confs of all mols, get the average runtime
       of all conformers for each molecule

    Parameters
    ----------
    sdfRef | str  | path+name of SDF file with times for all confs of all mols
    steps  | Bool | average number of steps instead of runtime seconds

    """
    
    # Open reference file.
    print("Opening SDF file %s" % sdfRef)
    ifsRef = oechem.oemolistream()
    ifsRef.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifsRef.open(sdfRef):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfRef)
    molsRef = ifsRef.GetOEMols()

    timeF = open(os.path.join(os.path.dirname(sdfRef),"timeAvgs.txt"), 'a')
    timeF.write("\nAnalyzing file: %s \n" % (sdfRef))
    if not steps: timeF.write("\nAverage runtime (sec) over all confs for each molecule\n")
    else: timeF.write("\nAverage number of steps over all confs for each molecule\n")

    # Grab all the times.
    for rmol in molsRef:
        if not steps:
            tmol = np.array(map(float, GetSDList(rmol, "runtime", method, basis)))
        else:
            tmol = np.array(map(float, GetSDList(rmol, "step", method, basis)))

        # exclude conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(tmol))
        for i in reversed(nanIndices): # loop in reverse to delete correctly
            tmol = np.delete(tmol, i)
        timeF.write( "%s\t%d confs\t\t%.3f\n" % (rmol.GetTitle(), tmol.size, np.mean(tmol)) )
    timeF.close()


def calcRelEne(sdfRef, method, basis, eFromOpt=False,outfn='relene.dat'):
    """

    WORDS

    Parameters--------UPDATE ME
    ----------
    sdf1 | str | path+name of SDF file with times for all confs of all mols
    sdf2 | str | should have same mols/confs as sdf1, likely diff coords/tags
    tag1 | str | data in this tag from sdf1 to be compared to sdf2
    tag2 | str | data in this tag from sdf2 to be compared to sdf1
    m1/b1| str | method/basis from sdf1. If m2/b2 is None, use same from m1/b1.
    verbose | bool | print information on each conformer

    For tags, see options in GetSDList function.

    """
    
    # Open file.
    print("Opening SDF file %s" % sdfRef)
    ifs1 = oechem.oemolistream()
    ifs1.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifs1.open(sdfRef):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sdf1)
    mols1 = ifs1.GetOEMols()

    # Determine SD tag from which to obtain energy.
    if eFromOpt: tagword = "QM opt energy"
    else: tagword = "QM spe"

    # Write description in output file.
    compF = open(os.path.join(os.path.dirname(sdfRef),outfn), 'w')
    compF.write("# Relative energies (kcal/mol) for file:\n# %s\n" % sdfRef) 
    compF.write("# using %s energies from %s/%s\n" % (tagword, method, basis))

    for imol in mols1:

        # Get absolute energies from the SD tags
        iabs = np.array(map(float, GetSDList(imol, tagword, method, basis)))

        # Get omega conformer number of first, for reference info
        # whole list can be used for matching purposes
        origids = GetSDList(imol, "original index")
        refconfid = origids[0]

        # For each, take relative energy to first then convert
        irel = iabs - iabs[0]
        irel = 627.5095*irel

        # Write output. This can be used for plotting, RMSD, etc.
        compF.write("\n# Mol title: %s\n# Energies relative to omega conf # %s\n"\
                   % (imol.GetTitle(), refconfid))
        for i in range(np.shape(irel)[0]):
            compF.write( "%s\t%.8f\n" % (origids[i], irel[i]) )
    compF.close()
    return True


### ------------------- Script -------------------

def main(**kwargs):
    if opt['time']:
        timeAvg(opt['filename'], opt['method'], opt['basisset'], steps=True)
        timeAvg(opt['filename'], opt['method'], opt['basisset'], steps=False)
    if opt['relene']:
        calcRelEne(opt['filename'], opt['method'], opt['basisset'], opt['efromopt'])


### ------------------- Parser -------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # analyzing results of an SDF file
    parser.add_argument("--time", action="store_true", default=False,
        help="If specified, get average number of steps and wall\
              clock time over all conformers of each molecule\
              reported in SDF file.")

    parser.add_argument("--relene", action="store_true", default=False,
        help="If specified, get lists of relative energies of all conformers\
 relative to first conformer in common that successfully completed its\
 QM calculation.")

    parser.add_argument("-f", "--filename",
        help="SDF file (with FULL path) to be processed.")
    parser.add_argument("-m", "--method",
        help="Name of QM method. Put this in 'quotes'.")
    parser.add_argument("-b", "--basisset",
        help="Name of QM basis set. Put this in 'quotes'.")
    parser.add_argument("--efromopt", action="store_true", default=False,
        help="If specified, get relative energies from optimization\
              instead of SPE (default).")

    args = parser.parse_args()
    opt = vars(args)

    ### Check file dependencies.

    if (opt['time'] or opt['relene']):
        if not os.path.exists(opt['filename']):
            raise parser.error("Input file %s does not exist." % opt['filename'])
        try: (opt['method'] and opt['basisset'])
        except NameError: raise parser.error("A method and basis set must be supplied for\
 either time or relative energy analysis.")

    main(**opt)

