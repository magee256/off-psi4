#!/usr/bin/env python

## This script takes a list of smiles as "filename.smi" and, for each smiles:
##  - conformers are generated with omega
##  - electrostatic clashes are resolved with Szybki minimization
##  - quick optimization is performed also with Szybki
##  - a new subdir is created in parent directory
## All conformers of all molecules are written out to "filename.sdf".

## Import and call smi2confs.smi2confs(arg1, arg2, arg3, arg4)

import os, sys
import openeye.oechem as oechem
import openeye.oeomega as oeomega
import openeye.oeszybki as oeszybki



### ------------------- Variables -------------------


#wdir = "/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-A"
#base = "AlkEthOH_chain_tiny"

### ------------------- Functions -------------------

def GenerateConfs(Mol):
    """
    Generate conformers of molecule from its SMILES string.

    Parameters
    ----------
    Mol:          OEChem molecule

    Returns
    -------
    molWithConfs: OEChem molecule with omega-generated conformers

    """
    molWithConfs = oechem.OEMol(Mol)
    omega = oeomega.OEOmega()
    maxConfs = 0
    omega.SetMaxConfs(maxConfs)
    omega.SetStrictStereo(False)
    omega.SetSampleHydrogens(True)
    omega.SetEnumNitrogen( oeomega.OENitrogenEnumeration_All)
    if not omega(molWithConfs):
        print( "omega failed on %s" % mol.GetTitle() )
        return None
    else:
        return molWithConfs

def ResolveBadClashes(Mol, Cfile):
    """
    Minimize conformers with severe steric interaction.

    Parameters
    ----------
    Mol:        single OEChem molecule (aka single conformer)
    Cfile:      string name of file to write output

    Returns
    -------
    boolean: True if completed successfully, False otherwise.

    """
    wfile = open(Cfile,'a')

    # set general energy options along with the single-point specification
    spSzybki = oeszybki.OESzybkiOptions()
    spSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    spSzybki.SetSolventModel(oeszybki.OESolventModel_Sheffield)
    spSzybki.SetRunType(oeszybki.OERunType_SinglePoint);
    # generate the szybki MMFF94 engine for single points
    szSP = oeszybki.OESzybki( spSzybki)
    # construct minimiz options from single-points options to get general optns
    optSzybki = oeszybki.OESzybkiOptions( spSzybki)
    # now reset the option for minimization
    optSzybki.SetRunType(oeszybki.OERunType_CartesiansOpt)
    # generate szybki MMFF94 engine for minimization
    szOpt = oeszybki.OESzybki( optSzybki)
    # add strong harmonic restraints to nonHs
    szOpt.SetHarmonicConstraints( 10.0)
    # construct a results object to contain the results of a szybki calculation
    szResults = oeszybki.OESzybkiResults()
    # work on a copy of the molecule
    tmpmol = oechem.OEMol( Mol)
    if not szSP(tmpmol, szResults):
        print( 'szybki run failed for %s' %  tmpmol.GetTitle() )
        return False
    Etotsp = szResults.GetTotalEnergy()
    Evdwsp = szResults.GetEnergyTerm(oeszybki.OEPotentialTerms_MMFFVdW)
    if Evdwsp > 35:
        if not szOpt(tmpmol, szResults):
            print( 'szybki run failed for %s' %  tmpmol.GetTitle() )
            return False
        Etot = szResults.GetTotalEnergy()
        Evdw = szResults.GetEnergyTerm(oeszybki.OEPotentialTerms_MMFFVdW)
        wfile.write( '%s resolved bad clash: initial vdW: %.4f ; '
                   'resolved EvdW: %.4f\n' % (tmpmol.GetTitle(),Evdwsp,Evdw) )
        Mol.SetCoords( tmpmol.GetCoords() )
    wfile.close()
    oechem.OESetSDData(Mol, oechem.OESDDataPair('MM Szybki Single Point Energy'\
, "%.12f" % szResults.GetTotalEnergy()))
    return True

def QuickOpt( Mol):
    """
    Fast MM optimization to whittle down number of conformers before QM.
    Default Szybki OEOptType type set to steepest descent (SD) based on
       preliminary comparisons.

    Parameters
    ----------
    Mol:        single OEChem molecule (aka single conformer)

    Returns
    -------
    boolean: True if completed successfully, False otherwise.

    """
    # set general energy options along with the run type specification
    optSzybki = oeszybki.OESzybkiOptions()
    optSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    optSzybki.SetSolventModel(oeszybki.OESolventModel_Sheffield)
    optSzybki.SetOptimizerType(oeszybki.OEOptType_SD)
    taglabel = 'MM Szybki SD Energy'

    # generate szybki MMFF94 engine for minimization
    szOpt = oeszybki.OESzybki( optSzybki)
    # construct a results object to contain the results of a szybki calculation 
    szResults = oeszybki.OESzybkiResults()
    # work on a copy of the molecule
    tmpmol = oechem.OEMol( Mol)
    if not szOpt(tmpmol, szResults):
        print( 'szybki run failed for %s' %  tmpmol.GetTitle() )
        return False
    Mol.SetCoords( tmpmol.GetCoords() )
    oechem.OESetSDData(Mol, oechem.OESDDataPair(taglabel, "%.12f" \
        % szResults.GetTotalEnergy()))
    return True


### ------------------- Script -------------------

def smi2confs(arg1, arg2, arg3=True, arg4=True):
    """
    From a file containing smiles strings, generate omega conformers,
       resolve steric clashes, do a quick MM opt, and write SDF output.

    Parameters
    ----------
    arg1: str - working directory containing .smi file
    arg2: str - name of the smiles file. E.g. "name.smi"
    arg3: boolean - Resolve steric clashes or not.
    arg4: boolean - QuickOpt or not.

    """
    wdir = arg1
    smiles = arg2
    sdfout = arg2.split('.')[0] + '.sdf'
    os.chdir(wdir)
    
    ### Read in smiles file.
    ifs = oechem.oemolistream()
    if not ifs.open(smiles):
        oechem.OEThrow.Warning("Unable to open %s for reading" % smiles)
    
    ### Open output file to write molecules.
    ofs = oechem.oemolostream()
    if os.path.exists(sdfout):
        #sys.exit("Output .sdf file already exists. Exiting.\n")
        print "Output .sdf file already exists. Exiting.\n"
        return
    if not ofs.open(sdfout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % sdfout)

    ### Output files detailing number of resolved clashes
    ###   and original number of conformers before MM opt.
    conffile = open('numOrigConfs.txt', 'w')
    
    ### For each molecule: label atoms, generate confs, resolve clashes, optimize.
    for smimol in ifs.GetOEMols():
        oechem.OETriposAtomNames(smimol)
        mol = GenerateConfs(smimol)
        conffile.write( "%s\t%s\n" % (mol.GetTitle(), mol.NumConfs()) )
    
        for i, conf in enumerate( mol.GetConfs()):
            print (mol.GetTitle(), i+1)
            ### Resolve bad clashes.
            if arg3:
                print "Resolving bad clashes..."
                if not ResolveBadClashes( conf, "numClashes.txt" ):
                    print('Resolving bad clashes failed for molecule %s \
conformer %d:' % (mol.GetTitle(),i+1) )
                    continue
            ### MM optimization.
            if arg4:
                print "Doing a quick MM (SD) optimization..."
                if not QuickOpt( conf):
                    print('Quick optimization failed for molecule %s \
conformer %d:' % (mol.GetTitle(),i+1) )
                    continue
        oechem.OEWriteConstMolecule(ofs, mol)
    
    ### Close files.
    ifs.close()
    ofs.close()
    conffile.close()

if __name__ == "__main__":
    smi2confs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
