#!/usr/bin/env python

# Usage: import, then call as 
#  timeAvg.timeAvg('/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-A/AlkEthOH_chain_tiny-optimized-minima.sdf')
import os
import openeye.oechem as oechem
import numpy as np


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
        if "JOB DID NOT FINISH" not in oechem.OEGetSDData(conf, "Note on opt."):
            SDList.append(oechem.OEGetSDData(conf, taglabel))
        else:
            SDList.append('nan')
    return SDList

### ------------------- Script -------------------


def timeAvg(sdfRef, steps=False):
    """

    For an SDF file with all confs of all mols, get the average runtime
       of all conformers for each molecule

    Parameters
    ----------
    sdfRef | str  | path+name of SDF file with times for all confs of all mols
    steps  | Bool | average number of steps instead of runtime seconds

    """
    
    # Open reference file.
    print("Opening SDF file ", sdfRef)
    ifsRef = oechem.oemolistream()
    ifsRef.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifsRef.open(sdfRef):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfRef)
    molsRef = ifsRef.GetOEMols()

    timeF = open(os.path.join(os.path.dirname(sdfRef),"timeAvgs.txt"), 'a')
    if not steps: timeF.write("\nAverage runtime (sec) over all confs for each molecule\n")
    else: timeF.write("\nAverage number of steps over all confs for each molecule\n")

    # Grab all the times.
    for rmol in molsRef:
        if not steps:
            tmol = np.array(map(float, GetSDList(rmol, "runtime", "mp2", "def2-sv(p)")))
        else:
            tmol = np.array(map(float, GetSDList(rmol, "step", "mp2", "def2-sv(p)")))
        timeF.write( "%s\t%s\n" % (rmol.GetTitle(), np.mean(tmol)) )
    timeF.close()

    print tmol


def compareSPEopt(sdf1, sdf2, tag1, tag2, m1, b1, m2=None, b2=None, verbose=False):
    """

    For an SDF file with all confs of all mols, get the average runtime
       of all conformers for each molecule

    Parameters
    ----------
    sdf1 | str | path+name of SDF file with times for all confs of all mols
    sdf2 | str | should have same mols/confs as sdf1, likely diff coords/tags
    tag1 | str | data in this tag from sdf1 to be compared to sdf2
    tag2 | str | data in this tag from sdf2 to be compared to sdf1
    m1/b1| str | method/basis from sdf1. If m2/b2 is None, use same from m1/b1.

    For tags, see options in GetSDList function.

    """
    
    # Open files.
    print("Opening SDF file ", sdf1)
    ifs1 = oechem.oemolistream()
    ifs1.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifs1.open(sdf1):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sdf1)
    mols1 = ifs1.GetOEMols()

    print("Opening SDF file ", sdf2)
    ifs2 = oechem.oemolistream()
    ifs2.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifs2.open(sdf2):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sdf1)
    mols2 = ifs2.GetOEMols()

    # Write description in output file.
    compF = open(os.path.join(os.path.dirname(sdf1),"comparison.txt"), 'a')
    compF.write("\nComparison 1 file: %s\nComparison 2 file: %s\n" % (sdf1, sdf2))

    # Figure out m2, b2
    if m2 is None: m2 = m1
    if b2 is None: b2 = b1

    for imol in mols1:
        jmol = mols2.next() # loop over both mols

        # Get absolute energies from the SD tags
        iabs = np.array(map(float, GetSDList(imol, tag1, m1, b1)))
        jabs = np.array(map(float, GetSDList(jmol, tag2, m2, b2)))

        # exclude conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(jabs))
        for i in reversed(nanIndices): # loop in reverse to delete correctly
            iabs = np.delete(iabs, i)
            jabs = np.delete(jabs, i)

        # For each, take relative energy to first 
        irel = iabs - iabs[0]
        jrel = jabs - jabs[0]

        # take RMSD of conformer energies for this particular mol
        dev = irel - jrel
        sqd = np.square(dev)
        mn = np.sum(sqd)/(np.shape(sqd)[0]-1)
        #print dev
        #print sqd
        #print mn
        rt = 627.5095*np.sqrt(mn)
        compF.write("\n%s\t%.4f\n" % (imol.GetTitle(), rt))
        if verbose:
            for i in range(np.shape(irel)[0]):
                compF.write( "%d\t%.8f\t%.8f\n" % (i, irel[i], jrel[i]) )
    compF.close()

