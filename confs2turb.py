#!/usr/bin/env python

## Usage: python confs2turb.py -i /path/and/filename.sdf

## Different from confs2psi.py since method, basisset, calc type, mem
##    should be specified in the templateOptions file of Turbomole setup.


import os
import argparse
import openeye.oechem as oechem
import subprocess as sp

### ------------------- Functions -------------------


def GetMolDetails(mol, label):
    """
    Parameters
    ----------
    mol: single OEChem conformer with coordinates
    label: string - name of the molecule. Can be an empty string.

    Returns
    -------
    optinfo: string - Turbomole input file for autoDefine.py
    xinfo:   string - XYZ format coordinates to feed into Turbomole's x2t

    """
    optinfo = ("$title %s" % label)
    optinfo += ("\n$charge %d" % oechem.OENetCharge( mol))
    optinfo += ("\n$end")

    xinfo = ("%d\n" % mol.NumAtoms())
    xyz = oechem.OEFloatArray(3)
    # get coordinates of each atom
    for atom in mol.GetAtoms():
        mol.GetCoords( atom, xyz)
        xinfo+=( '\n  %s %10.4f %10.4f  %10.4f' \
                  %(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()),
                  xyz[0], xyz[1], xyz[2]) )
    return optinfo, xinfo


### ------------------- Script -------------------

def confs2turb(insdf):
    """
    Parameters
    ----------
    insdf:  string - PATH+name of SDF file

    """
    homedir = os.getcwd()
    p = sp.call('module load turbomole/7.1/intel', shell=True)
    
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
            # change into subdirectory to use x2t
            subdir = os.path.join(homedir,"%s/%s" % (mol.GetTitle(), i+1))
            if not os.path.isdir(subdir):
                os.makedirs(subdir)
            os.chdir(subdir)

            # write out relevant files 
            label = mol.GetTitle()+'_'+str(i+1)
            ofile = open('options','w')
            xfile = open('input.xyz','w')
            optinfo, xinfo = GetMolDetails(conf,label)
            ofile.write(optinfo)
            xfile.write(xinfo)
            ofile.close()
            xfile.close()

            # run x2t
            p=sp.Popen('x2t input.xyz > coord',shell=True)
            p.wait()
            #os.chdir(wdir) # i don't think i need this?
            
    ifs.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script generates, for \
 each conformer, a Turbomole-style coord file as well as an options\
 file for use with autoDefine.py which automates define process of Turbomole.\
 Options file contains title and charge of mol.\
 x2t is run for each coord file to generate Turbomole coordinates.')

    parser.add_argument('-i','--infile', help='Input file with mols and confs\
 for each mol. Include the full path with filename.')

    args = parser.parse_args()
    confs2turb(args.infile)
