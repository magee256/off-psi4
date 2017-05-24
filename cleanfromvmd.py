#!/usr/bin/env python

# Usage: python sdfMol2.py path/with/refFile path/with/filteredFile path/with/newFiltfile

import os,sys
import openeye.oechem as oechem

def convertSDFfile(reffile, filtfile, writeout):
    refifs = oechem.oemolistream()
    filtifs = oechem.oemolistream()
    ofs = oechem.oemolostream()
    
    ### Read in reference file, but don't need its old conformers
    if not refifs.open(reffile):
        oechem.OEThrow.Warning("Unable to open %s for reading" % reffile)
        return

    ### Read in filtered file and distinguish each molecule's conformers
    filtifs.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not filtifs.open(filtfile):
        oechem.OEThrow.Warning("Unable to open %s for reading" % filtfile)
        return

    ### Open outstream file.
    if os.path.exists(writeout):
        print("File already exists: %s. Skip getting results.\n" % (finsdf))
        return
    if not ofs.open(writeout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % writeout)

    ### Loop and write molecules. (though refifs should only have ONE mol)
    for rmol in refifs.GetOEMols():
        for fmol in filtifs.GetOEMols():
            for i, conf in enumerate( fmol.GetConfs()):
                rmol.SetCoords(conf.GetCoords())
                oechem.OEWriteConstMolecule(ofs, rmol)

    refifs.close()
    filtifs.close()
    ofs.close()

if __name__ == "__main__":
    convertSDFfile(sys.argv[1], sys.argv[2], sys.argv[3])
