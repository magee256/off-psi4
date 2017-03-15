#!/usr/bin/env python

## This script generates Psi4 input files for each conf of each molecule.
## Specifying SPE=True will write inputs for single point energy calcns,
## else the default input file type is for geometry optimization.

## Import and call confs2psi.confs2psi(insdf, method, basis)


import os, sys
import openeye.oechem as oechem
import shutil


### ------------------- Functions -------------------


def MakePSI4Input(mol, label, method, basisset, SPE=False, mem=None):
    """
    Parameters
    ----------
    mol: single OEChem conformer with coordinates
    label: string - name of the molecule. Can be an empty string.
    method: string - specification of method (see Psi4 website for options)
    basisset: string - specification of basis set
    SPE: boolean - False (default) for geom opt. True for single point E calcns
    mem: string - specify Psi4 job memory. E.g. "2 Gb" "2000 Mb" "2000000 Kb"

    Returns
    -------
    inputstring: string - containing contents of whole input file for this conf

    """
    inputstring = ""
    xyz = oechem.OEFloatArray(3)
    # specify memory requirements, if defined
    if mem != None:
        inputstring += "memory %s\n" % mem
    inputstring+=( 'molecule %s {\n' % label )
    # charge and multiplicity; multiplicity hardwired to singlet (usually is)
    netCharge = oechem.OENetCharge( mol)
    inputstring+=( '  %s 1' % netCharge )
    # get coordinates of each atom
    for atom in mol.GetAtoms():
        mol.GetCoords( atom, xyz)
        inputstring+=( '\n  %s %10.4f %10.4f  %10.4f' \
                       %(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()),
                       xyz[0], xyz[1], xyz[2]) )
    inputstring+=( '\n  units angstrom')
    # explicitly specify MP2 RI-auxiliary basis for Ahlrichs basis set
    # http://www.psicode.org/psi4manual/master/basissets_byfamily.html
    if method.lower()=='mp2' and 'def' in basisset and basisset.lower()!='def2-qzvpd':
        inputstring+=('\n}\n\nset basis %s' % (basisset))
        inputstring+=('\nset df_basis_mp2 %s-ri' % (basisset))
        inputstring+=('\nset freeze_core True')
    else:
        inputstring+=('\n}\n\nset basis %s' % (basisset))
        inputstring+=('\nset freeze_core True')
    # specify command for type of calculation
    if SPE is False:
        inputstring+=('\noptimize(\'%s\')' % (method))
    else:
        inputstring+=('\nenergy(\'%s\')' % (method))
    return inputstring


### ------------------- Script -------------------

def confs2psi(insdf, method, basis, spe=False, memory=None):
    """
    Parameters
    ----------
    insdf:  string - PATH+name of SDF file
    method: string - method. E.g. "mp2"
    basis:  string - basis set. E.g. "def2-sv(p)"
    spe:    boolean. True for single point energy calcns, False for geom opt.
            default option is False.
    memory: string - memory specification. Psi4 default is 256 Mb. E.g. "1.5 Gb"

    """
    wdir, fname = os.path.split(insdf)
    os.chdir(wdir)
    
    ### Read in .sdf file and distinguish each molecule's conformers
    ifs = oechem.oemolistream()
    ifs.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifs.open(insdf):
        oechem.OEThrow.Warning("Unable to open %s for reading" % insdf)
        return
    
    ### For each molecule: for each conf, generate input
    for mol in ifs.GetOEMols():
        print(mol.GetTitle(), mol.NumConfs())
        for i, conf in enumerate( mol.GetConfs()):
            # change into subdirectory ./mol/conf/
            subdir = os.path.join(wdir,"%s/%s" % (mol.GetTitle(), i+1))
            if not os.path.isdir(subdir):
                os.makedirs(subdir)
            if os.path.exists(os.path.join(subdir,'input.dat')):
                print "Input file (\"input.dat\") already exists. Skipping.\n"
                continue
            label = mol.GetTitle()+'_'+str(i+1)
            ofile = open(os.path.join(subdir,'input.dat'), 'w')
            ofile.write(MakePSI4Input( conf, label, method, basis, spe, memory))
            ofile.close()
    ifs.close()

def prep(source, sink):

    """
               untested 1-17-17
    Copies source file to be sink file. Check for already existing sink file.
    Makes directory(ies) for sink file if necessary.

    Parameters
    ---------- 
    source: string, filename with absolute path
    sink:   string, filename with absolute path

    """

    newdir = os.path.dirname(sink)

    # check for existing file
    if os.path.exists(sink):
        print("File already exists: %s. Please provide a new filename.\n" % (sink))
        return

    # make new directory if needed
    if not os.path.isdir(newdir):
        os.makedirs(newdir)

    # copy file from olddir to newdir
    try:
        shutil.copy2(source, sink)
    except IOError, e:
        print(e)
        return
