
## TODO: make so that cc-pVTZ == cc-pvtz == whatever else
## TODO: plot from verbose file, read in mol title, plot like with other case, large dots

# Note: If you see error: "ValueError: could not convert string to float:"
#   check to make sure that the value for SDF tag is correct.
import os
import openeye.oechem as oechem
import numpy as np
import argparse

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

def calcRelEne(sdf1, sdf2, m1, b1, m2, b2, spe1=True, spe2=True, verbose=True,outfn='relene-rmsd.dat'):
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

    def prelim(sdfRef,spe):
        # Open file.
        print("Opening SDF file %s" % sdfRef)
        ifs1 = oechem.oemolistream()
        ifs1.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifs1.open(sdfRef):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % sdf1)
        mols = ifs1.GetOEMols()
    
        # Determine SD tag from which to obtain energy.
        if spe: tagword = "QM spe"
        else: tagword = "QM opt energy"
        return mols, tagword

    mols1, tag1 = prelim(sdf1,spe1)
    mols2, tag2 = prelim(sdf2,spe2)

    # Write description in output file.
    compF = open(os.path.join(os.path.dirname(sdf1),outfn), 'w')
    compF.write("# RMSD of relative energies (kcal/mol):\n")
    compF.write("# File1: %s\n" % sdf1) 
    compF.write("#   using %s energies from %s/%s\n" % (tag1, m1, b1)) 
    compF.write("# File2: %s\n" % sdf2) 
    compF.write("#   using %s energies from %s/%s\n" % (tag2, m2, b2)) 

    for imol in mols1:
        jmol = mols2.next()

        # Get absolute energies from the SD tags
        print GetSDList(imol, tag1, m1, b1)
        iabs = np.array(map(float, GetSDList(imol, tag1, m1, b1)))
        jabs = np.array(map(float, GetSDList(jmol, tag2, m2, b2)))

        # Get omega conformer number of first, for reference info
        # whole list can be used for matching purposes
        originum = GetSDList(imol, "original index")
        origjnum = GetSDList(jmol, "original index")


        # exclude conformers for which job did not finish (nan)
        # check for file1
        nanIndices = np.argwhere(np.isnan(iabs))
        iabs = np.delete(iabs,nanIndices)
        jabs = np.delete(jabs,nanIndices)
        originum = np.delete(np.asarray(originum),nanIndices)
        origjnum = np.delete(np.asarray(origjnum),nanIndices)
        # same check for file2
        nanIndices = np.argwhere(np.isnan(jabs))
        iabs = np.delete(iabs,nanIndices)
        jabs = np.delete(jabs,nanIndices)
        originum = np.delete(np.asarray(originum),nanIndices)
        origjnum = np.delete(np.asarray(origjnum),nanIndices)

        # Take relative energy to first conf
        irel = iabs - iabs[0]
        jrel = jabs - jabs[0]


        # take RMSD of conformer energies for this particular mol
        dev = irel - jrel
        sqd = np.square(dev)
        mn = np.sum(sqd)/(np.shape(sqd)[0]-1)
        rt = 627.5095*np.sqrt(mn)

        # Write output. This can be used for plotting, RMSD, etc.
        compF.write("\n# Mol %s, RMSD = %.5f kcal/mol\n" % (imol.GetTitle(), rt))
        compF.write("# Energies relative to omega conf #%s\n"\
                   % (originum[0]))
        if verbose:
            # convert relative energies from Hartrees to kcal/mol
            irel = 627.5095*irel
            jrel = 627.5095*jrel
            for i in range(np.shape(irel)[0]):
                compF.write( "%s\t%.8f\t%.8f\n" % (originum[i], irel[i], jrel[i]) )
    compF.close()
    return True

sdf1='/data12/cmf/limvt/qm_AlkEthOH/pipeline/comparison-ab/verse/divrefine-221-spe1.sdf'
sdf2='/data12/cmf/limvt/qm_AlkEthOH/pipeline/comparison-ab/verse/divrefine-221-opt2.sdf'
#sdf2='/data12/cmf/limvt/qm_AlkEthOH/pipeline/comparison-ab/verse/divrefine-221-spe2.sdf'
#sdf2='/data12/cmf/limvt/qm_AlkEthOH/pipeline/comparison-ab/verse/divrefine-221-spe3.sdf'
m1 = 'b3lyp-d3mbj'
b1 = 'def2-tzvp'
m2 = 'b3lyp-d3mbj'
b2 = 'def2-tzvp'
#m2 = 'mp2'
#b2 = 'cc-pVTZ'
#m2 = 'PBE0'
#b2 = '6-311G**'

calcRelEne(sdf1, sdf2, m1, b1, m2, b2, spe1=True, spe2=False,outfn='relene-rmsd_10.dat')

# dictionary
# -  (1) filename, (2) method, (3) basis, (4) boolean if spe, (5) boolean if is reference 
# how do dictionaries work
