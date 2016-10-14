#!/usr/bin/env python

## This script takes an SDF file and, for each molecule:
##  - conformers are compared by energy
##  - conformers are compared by RMSD
## in order to roughly filter out duplicate minima and keep unique ones.
## Filtered conformers for all molecules are written out in SDF file.

## Import and call filterConfs.filterConfs(arg1, arg2, arg3)

import re
import os, sys, glob
import openeye.oechem as oechem


### ------------------- Variables -------------------

### Parameters for distinguishing minima & OERMSD
automorph = True
heavyOnly = False
overlay = True
thresE = 5.E-4 # declare confs diff & skip RMSD comparison above this threshold
thresRMSD = 0.2 # above this threshold (Angstrom), confs are "diff" minima


### ------------------- Functions -------------------

def IdentifyMinima(Mol,Taglabel,ThresholdE,ThresholdRMSD):
    """
    For a molecule's set of conformers computed with some level of theory, 
        whittle down unique conformers based on energy and RMSD. 

    Parameters
    ----------
    Mol           OEChem molecule with all of its conformers
    Taglabel      string name of the SD tag in this molecule
    ThresholdE    float value for abs(E1-E2), below which 2 confs are "same"
        Units are hartrees (default output units of Psi4)
    ThresholdR    float value for RMSD, below which 2 confs are "same"
        Units are in Angstrom (Psi4 default)

    Returns
    -------
    boolean True if successful filter + delete. False if there's only 
        one conf and it didn't optimize, or something else funky.

    """
    confsToDel = set() # declare an empty set (unordered) for confs to delete

    # if there's only 1 conf, has SDData ==> True, not ==> False
    if Mol.NumConfs()==1:
        testmol = next(Mol.GetConfs())
        if not oechem.OEHasSDData(testmol, Taglabel):
            return False
        else:
            return True

    # Loop over conformers twice (NxN diagonal comparison of RMSDs)
    for confRef in Mol.GetConfs():
        print(" ~ Reference: %s conformer %d" % (Mol.GetTitle(), confRef.GetIdx()+1))
        # delete cases that don't have energy (opt not converged; or other)
        if not oechem.OEHasSDData(confRef, Taglabel):
            confsToDel.add(confRef.GetIdx())
            continue
        refE = float(oechem.OEGetSDData(confRef,Taglabel))
    
        for confTest in Mol.GetConfs():
            # upper right triangle comparison
            if confTest.GetIdx() <= confRef.GetIdx(): 
                continue 
            # skip cases already set for removal
            if confTest.GetIdx() in confsToDel:
                continue
            # delete cases that don't have energy
            if not oechem.OEHasSDData(confTest, Taglabel):
                confsToDel.add(confTest.GetIdx())
                continue

            testE = float(oechem.OEGetSDData(confTest,Taglabel))
            # if MM (not Psi4) energies, convert absERel to Hartrees
            if 'mm' in Taglabel.lower():
                absERel = abs(refE-testE)/627.5095
            else: 
                absERel = abs(refE-testE)
            # if energies are much diff., confs are diff, skip.
            if absERel > ThresholdE:
                #print "Skip from energy difference ",absERel
                continue
            # for the confs with similar E, see if they are diff with RMSD
            rmsd = oechem.OERMSD(confRef,confTest,automorph,heavyOnly,overlay)
            # if RMSD < thresholdRMSD, must be same conf, tag to delete.
            if rmsd < ThresholdRMSD:
                confsToDel.add(confTest.GetIdx()) 
            
    # for the same molecule, delete tagged conformers
    print Mol.GetTitle(), " original number of conformers: ",Mol.NumConfs()
    for conf in Mol.GetConfs():
        if conf.GetIdx() in confsToDel:
            print 'Removing %s conformer index %d' \
                  % (Mol.GetTitle(),conf.GetIdx())
            if not Mol.DeleteConf(conf):
                oechem.OEThrow.Fatal("Unable to delete %s GetIdx() %d" \
                                  % (Mol.GetTitle(), conf.GetIdx()))
    return True


### ------------------- Script -------------------

def filterConfs(arg1, arg2, arg3):
    wdir = arg1
    rmsdfile = arg2
    tag = arg3

    os.chdir(wdir)
    numConfsF = open(os.path.join(wdir,"numFiltConfs.txt"), 'a')
    numConfsF.write(tag+"\n")

    # Open file to be processed.
    rmsd_ifs = oechem.oemolistream()
    if not rmsd_ifs.open(rmsdfile):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % rmsdfile)
    rmsd_ifs.SetConfTest( oechem.OEAbsoluteConfTest() )
    rmsd_molecules = rmsd_ifs.GetOEMols()
    
    # Open outstream file.
    rmsdout = arg2.split('.')[0] + "-minima.sdf"
    rmsd_ofs = oechem.oemolostream()
    if os.path.exists(rmsdout):
        print("%s output file already exists. Skip filtering.\n" % rmsdout)
        return
    if not rmsd_ofs.open(rmsdout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % rmsdout)

    # Identify minima and write output file.         
    for mol in rmsd_molecules:
        if IdentifyMinima(mol,tag,thresE,thresRMSD):
            numConfsF.write( "%s\t%s\n" % (mol.GetTitle(), mol.NumConfs()) )
            oechem.OEWriteConstMolecule(rmsd_ofs, mol)
        else:
            numConfsF.write( "%s\t0\n" % (mol.GetTitle()))
    rmsd_ifs.close()
    numConfsF.close()
    rmsd_ofs.close()

    print("Done filtering %s.\n" % (rmsdfile))
# done

