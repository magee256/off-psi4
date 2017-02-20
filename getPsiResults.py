#!/usr/bin/env python

import re
import os, sys, glob
import openeye.oechem as oechem



### ------------------- Functions -------------------

def GetTime(filename):
    """ 
    Get wall-clock time from Psi4 file. If multiple times are present,
        the average will be taken. Used in CompareTimes(...) function.

    Parameters
    ----------
    filename: string name of the timefile. E.g. "timer.dat"

    Returns
    -------
    time: float of the average wall-clock time of a single timefile

    """
    with open(filename) as fname:
        times = []
        for line in fname:
            if "Wall Time:" in line:
                times.append(float(line.split()[2]))
    time = sum(times) / float(len(times))
    return time

def ProcessOutput(filename, Props, spe=False):
    """ 
    Go through output file and get level of theory (method and basis set),
        number of optimization steps, initial and final energies, and
        optimized coordinates. Returns all this information in a dictionary
        that was passed to this function. 

    Parameters
    ----------
    filename: string name of the output file. E.g. "output.dat"
    Props: dictionary where all the data will go. Can be empty or not.
    spe: Boolean - are the Psi4 results of a single point energy calcn?

    Returns
    -------
    Props: dictionary with summarized data from output file.
           keys are: basis, method, numSteps, initEnergy, finalEnergy, coords

    """
    try:
        f = open(filename,"r")
    except IOError:
        sys.exit("No %s file found in directory of %s" \
                 % ( filename,os.getcwd() ))


    rough = []
    coords = []
    lines = f.readlines()
    it = iter(lines)

    if spe: # Process single point energy calcns differently
        for line in it:
            if "set basis" in line:
                Props['basis'] = line.split()[2]
            if "energy(" in line:
                Props['method'] = line.split('\'')[1]
            if "Total Energy =" in line:
                Props['finalEnergy'] = float(line.split()[3])
        return Props

    # Loop through file to get method, basis set, numSteps, energies, coords
    for line in it:
        if "set basis" in line:
            Props['basis'] = line.split()[2]
        if "optimize(" in line:
            Props['method'] = line.split('\'')[1]
        if "Optimization is complete" in line:
            Props['numSteps'] = line.strip().split(' ')[5]
            for _ in xrange(8):
                line = next(it)
            Props['initEnergy'] = float(line.split()[1])
        if "Final energy" in line:
            Props['finalEnergy'] = float(line.split()[3])
            line = next(it) # "Final (previous) structure:"
            line = next(it) # "Cartesian Geometry (in Angstrom)"
            line = next(it) # Start of optimized geometry
            while "Saving final" not in line:
                rough.append(line.split()[1:4])
                line = next(it)

    # Convert the 2D (3xN) coordinates to 1D list of length 3N (N atoms).
    for atomi in rough:
        coords += [float(i) for i in atomi]
    Props['coords'] = coords
    f.close()
    return Props

    
def SetOptSDTags(Conf, Props, spe=False):
    """
    WORDS WORDS WORDS

    Parameters
    ----------
    Conf:       Single conformer from OEChem molecule
    Props:      Dictionary output from ProcessOutput function.
                Should contain the keys: basis, method, numSteps,
                initEnergy, finalEnergy, coords, time
    spe:        Boolean - are the Psi4 results of a single point energy calcn?

    """


    # get level of theory for setting SD tags
    method = Props['method']
    basisset = Props['basis']
    
    # check that finalEnergy is there. if not, opt probably did not finish
    # make a note of that in SD tag
    if not 'finalEnergy' in Props:
        if not spe: oechem.OEAddSDData(Conf, "Note on opt.", \
 "JOB DID NOT FINISH (opt %s/%s)" % (method, basisset))
        else: oechem.OEAddSDData(Conf, "Note on opt.",\
 "JOB DID NOT FINISH (opt %s/%s)" % (method, basisset))
        return

    # Set new SD tag for conformer's final energy
    if not spe: taglabel = "QM Psi4 Final Opt. Energy (Har) %s/%s" % (method, basisset)
    else: taglabel = "QM Psi4 Single Pt. Energy (Har) %s/%s" % (method, basisset)
    oechem.OEAddSDData(Conf, taglabel, str(Props['finalEnergy']))

    # Set new SD tag for wall-clock time
    if not spe: taglabel = "QM Psi4 Opt. Runtime (sec) %s/%s" % (method, basisset)
    else: taglabel = "QM Psi4 Single Pt. Runtime (sec) %s/%s" % (method, basisset)
    oechem.OEAddSDData(Conf, taglabel, str(Props['time']))

    if spe: return # stop here if SPE

    # Set new SD tag for original conformer number
    # !! Opt2 files should ALREADY have this !! Opt2 index is NOT orig index!
    taglabel = "Original omega conformer number"
    if not oechem.OEHasSDData(Conf, taglabel): # only add tag if not existing
        oechem.OEAddSDData(Conf, taglabel, str(Conf.GetIdx()+1))

    # Set new SD tag for numSteps of geom. opt.
    taglabel = "QM Psi4 Opt. Steps %s/%s" % (method, basisset)
    oechem.OEAddSDData(Conf, taglabel, str(Props['numSteps']))

    # Set new SD tag for conformer's initial energy
    taglabel = "QM Psi4 Initial Opt. Energy (Har) %s/%s" % (method, basisset)
    oechem.OEAddSDData(Conf, taglabel, str(Props['initEnergy']))


### ------------------- Script -------------------

def getPsiResults(wdir, origsdf, finsdf, spe=False, timefile=None, psiout=None):
    """
    Parameters
    ----------
    wdir: string - directory containing (1) all confs' jobs, (2) orig sdf file,
                   and the (3) soon-to-be-generated final output sdf file.
    origsdf:  string - full name of original pre-opt SDF file. E.g. "name.sdf"
    finsdf:   string - full name of final SDF file with optimized results.
    spe:     Boolean - are the Psi4 results of a single point energy calcn?
    timefile: string - name of the Psi4 timer files. Default is "timer.dat"
    psiout:   string - name of the Psi4 output files. Default is "output.dat"

    """
    if timefile is None:
        timefile = "timer.dat"
    if psiout is None:
        psiout = "output.dat"

    os.chdir(wdir)
    # Read in .sdf file and distinguish each molecule's conformers
    os.chdir(wdir)
    ifs = oechem.oemolistream()
    ifs.SetConfTest( oechem.OEAbsoluteConfTest() )
    if not ifs.open(origsdf):
        oechem.OEThrow.Warning("Unable to open %s for reading" % origsdf)
        quit()
    molecules = ifs.GetOEMols()
    
    # Open outstream file.
    writeout = os.path.join(wdir,finsdf)
    write_ofs = oechem.oemolostream()
    if os.path.exists(writeout):
        print("File already exists: %s. Skip getting results.\n" % (finsdf))
        return
    if not write_ofs.open(writeout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % writeout)
    
    # For each conformer, process output file and write new data to SDF file
    for mol in molecules:
        print("===== %s =====" % (mol.GetTitle()))
        for j, conf in enumerate( mol.GetConfs()):
            props = {} # dictionary of data for this conformer 
            # change into subdirectory ./mol/conf/
            subdir = os.path.join(wdir,"%s/%s" % (mol.GetTitle(), j+1))
            if not os.path.isdir(subdir):
                print("*** No directory found for %s ***" % (subdir))
                continue
            os.chdir(subdir)
            # Get wall clock time of the job
            try:
                props['time'] = GetTime(timefile)
            except IOError:
                props['time'] = "timer.dat file not found"
                pass
            # process output and get dictionary results
            #if not spe: props = ProcessOutput(psiout, props)
            #else: props = ProcessOutput(psiout, props, spe=True)
            props = ProcessOutput(psiout, props, spe)
            # Set last coordinates from optimization
            #    skip SetCoords if coords missing
            if 'coords' in props and len(props['coords']) != 0 :
                conf.SetCoords(oechem.OEFloatArray(props['coords']))
            # Set SD tags for this molecule
            #if not spe: SetOptSDTags(conf, props)
            #else: SetOptSDTags(conf, props, spe)
            SetOptSDTags(conf, props, spe)
            # Write output file
            oechem.OEWriteConstMolecule(write_ofs, conf)
    ifs.close()
    write_ofs.close()
    try:
        return props['method'], props['basis']
    except KeyError:
        return None, None 

if __name__ == "__main__":
    getPsiResults(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

