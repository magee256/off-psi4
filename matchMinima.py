#!/usr/bin/env python

import re
import os, sys, glob
import pickle
import collections
import openeye.oechem as oechem
import numpy as np
from numpy import nan
from itertools import product # letter labels for minima
import matplotlib.pyplot as plt
import matplotlib as mpl
import shutil

### ------------------- Variables -------------------

### Parameters for connecting minima by RMSD
automorph = True
heavyOnly = False
overlay = True

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
        taglabel = "QM Psi4 Final Opt. Energy %s/%s" % (Method, Basisset)

    if Property=="MM opt energy":
        taglabel = "MM Szybki Newton Energy"

    if Property=="original index":
        taglabel = "Original omega conformer number"

    SDList = []
    for j, conf in enumerate( Mol.GetConfs() ):
        SDList.append(oechem.OEGetSDData(conf, taglabel))
    return SDList

def Compare2Mols(Rmol, Qmol):
    # 1D list of length N = number of minimas in reference file.
    # Index of qmol that best matches with each rmol minima.
    MolIndices = []
    rmsdConf = [] # for this conf, get indices of matching minima of other theories
    for Rconf in Rmol.GetConfs():
        print ">>>> Matching conformers to minima: %d <<<<"\
            % (Rconf.GetIdx()+1)
        rsublist = [] # list of RMSDs for this level of theory
        for Qconf in Qmol.GetConfs():
            rms = oechem.OERMSD(Rconf,Qconf,automorph,heavyOnly,overlay)
            #print Qconf.GetIdx()+1, rms
            rsublist.append(rms)

        # get the index for the minimum value RMSD
        thisMin=[i for i, j in enumerate(rsublist) if j == min(rsublist)][0]
        # deem it a match if the minimum RMSD is <= 0.1 A
        if rsublist[thisMin] <= 0.1:
            rmsdConf.append(thisMin)
        else:
            rmsdConf.append(None)
    MolIndices.append(rmsdConf)
    return MolIndices

def PlotMolMinima(MolName, RefNumConfs, MinimaE, Xticklabels, RefFile, Wdir, Stag=False):
    # stagger plots to better see line overlap. works best with few (<5?) numFiles.
    # assumes minimaE is generated correctly and all sublists in minimaE overlap.

    # minimaE for ONE molecule

    numFiles = len(MinimaE)
    ### Flatten this 2D list into a 1D to find min and max for plot
    flatten = [item for sublist in MinimaE for item in sublist]
    floor = min(flatten)
    ceiling = max(flatten)
    ystep = (ceiling - floor)/9
    ystep = round(ystep * 2) / 2 # round the step to nearest 0.5

    ### Stagger each of the component files of minimaE for ease of viewing.
    if Stag==True:
        tempMinimaE = []
        for i, fileE in enumerate(MinimaE):
            tempMinimaE.append([x+i/2. for x in fileE])
        MinimaE = tempMinimaE
        ceiling = ceiling + numFiles

    ### Figure labels.
    plttitle="Relative Energies of %s Minima" % (MolName)
    plttitle+="\nby Reference File %s" % (RefFile)
    ylabel="Relative energy (kcal/mol)"
    figname = os.path.join(Wdir,"%s-minimaE.png" % (MolName))
    
    ### Set x-axis values and labels. Can either be letter or number, not both.
    # LETTER LABELS
    letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ' # label x-axis by letter
    rpt = (len(MinimaE[0])/26)+1
    xlabs =[''.join(i) for i in product(letters, repeat=rpt)][:RefNumConfs]
    # NUMBER LABELS
    xlabs = range(len(MinimaE[0]))
    
    ### Set this up for grid.
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xticks(np.arange(-1,RefNumConfs+1,2))
    
    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle,fontsize=16)
    plt.ylabel(ylabel,fontsize=14)
    plt.xlabel("minimum",fontsize=14)
    plt.xticks(range(RefNumConfs),xlabs,fontsize=12)
    plt.yticks(fontsize=12)
    
    ### Plot the data.
    colors = mpl.cm.rainbow(np.linspace(0, 1, numFiles))
    markers = ["x","^","8","d","o","s","*","p","v","<",".","+",">","D"]
    colors = ['r','y','b']
    markers = [".","^","d"]
    #for i in reversed(range(numFiles)):
    for i, FileE in enumerate(MinimaE):
        xi = range(RefNumConfs)
        #yi = [item[i] for item in MinimaE]
        yi = FileE
        plt.plot(xi,yi,color=colors[i],label=Xticklabels[i],\
            marker=markers[i],markersize=9)
    
    ### Add legend and set plot limits.
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.xlim(-1,RefNumConfs+1)
    # y axis limits: min, max, step
    ax.set_yticks(np.arange(floor-2,ceiling+2,ystep))
    plt.grid()
   
    plt.savefig(figname,bbox_inches='tight')
    plt.show()
### ------------------- Script -------------------

def prep(arg1, arg2, arg3, arg4=None):
    """
    Combine SDF files with optimized geometries/energies into
       directory for minima matching. 

    Parameters
    ----------
    arg1: str - directory containing SDF file to be analyzed
    arg2: str - name of SDF file to be analyzed
    arg3: str - matchMinima directory where all SDF files will be analyzed
    arg4: str - optional parameter for SDF file to be renamed in destination

    """
    origdir = arg1
    source = arg2
    destdir = arg3
    if arg4 != None:
        sink = arg4
    else:
        sink = arg1.split('/')[-1] + '.sdf'
    # check or make destination directory
    if not os.path.isdir(destdir):
        os.makedirs(destdir)
    # check if sink file already exists in destination
    if os.path.exists(os.path.join(destdir,sink)):
        print "Match file already exists. Skipping. %s" %\
           os.path.join(destdir,sink)
        return

    try:
        shutil.copy2(os.path.join(origdir, source), os.path.join(destdir, sink))
    except IOError, e:
        print(e)
        return
    print "matchMinima file prep complete for %s." %\
        os.path.join(destdir, sink)

def matchMinima(wdir, sdfRef, sdfList):
    """
    For list of SDF files, match the conformer minima to those of the reference
       SDF file. Ex. Conf #1 of reference matches with conf #7 of file1 and
       with conf #18 of file2.
    Each SDF query file should have the same molecules in the same order
       directory for minima matching. Number of confs per mol can differ.

    Parameters
    ----------
    wdir: str - working directory containing SDF files to be analyzed
    sdfRef: str - name of the SDF file which is to be used as reference for all mols
    sdfList: str list - list of the SDF file names to be analyzed.
          This list should include the reference SDF file (sdfRef).

    """
    
    numFiles = len(sdfList)
    
    os.chdir(wdir)
    numMols = 0
    allIndices = [] # for M mols, N reference minima of each mol, P matching indices for each ref minimia
    elists = [] # 2D list: K mols per file x J numFiles
    refNumConfs = [] # number of conformers for each mol in reference file
    molNames = [] # name of each molecule. for plotting.
    
    for i, sdfQuery in enumerate(sdfList):
        # Open reference file.
        print "Opening reference file ", sdfRef
        ifsRef = oechem.oemolistream()
        ifsRef.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifsRef.open(sdfRef):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfRef)
        molsRef = ifsRef.GetOEMols()
    
        # Open query file.
        print "Opening query file ", sdfQuery
        ifsQuery = oechem.oemolistream()
        ifsQuery.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifsQuery.open(sdfQuery):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfQuery)
        molsQuery = ifsQuery.GetOEMols()

        for rmol in molsRef:
            qmol = molsQuery.next() # loop over both mols.
            molNames.append(rmol.GetTitle())
    
            # get energies for plotting relative energies
            elists.append(GetSDList(qmol, "QM opt energy", "mp2", "def2-sv(p)"))
    
            # Skip minmatch if this is the same as the reference theory.
            # Do this here because need to fill elists, get refNumConfs int, 
            #   add placeholder list for molIndices.
            if(sdfQuery == sdfRef):
                print "\n\nSkipping comparison against self."
                refNumConfs.append(rmol.NumConfs()) # num of confs in ref
                allIndices.append([-1]*rmol.NumConfs())
                continue
            molIndices = Compare2Mols(rmol, qmol)[0] # [[double-listed]] ####
            allIndices.append(molIndices) ####
            numMols = numMols + 1
    numMols = numMols / (numFiles-1) # avoid repeat counting
    molNames = molNames[:numMols]

#    print "\n\nmolNames\n",molNames
#    print "\n\nrefNumConfs\n",refNumConfs
#    print "\n\nallIndices\n",allIndices
#    print "\n\nelists\n",elists


    pickle.dump([molNames, refNumConfs,allIndices,elists], open(os.path.join(wdir,'matchMinima.pickle'), 'wb'))
#    molNames, refNumConfs, allIndices, elists = pickle.load(open(os.path.join(wdir,'matchMinima.pickle'), 'rb'))

    ### List is now: [[file1 mol1] ... [file1 molN] [file2 mol1] ... [file2 molN]]
    ### Reorganize to [[file1 mol1] [file2 mol1] ... [file1 molN] [file2 molN]]
    if len(allIndices) != (numFiles * numMols):
        oechem.OEThrow.Fatal("Length of allIndices list is %d. It should\
be (numFiles * numMols = %d) Exiting." % (len(allIndices), numFiles * numMols))
    allIndices = [allIndices[i::numMols] for i in range(numMols)] 
    elists = [elists[i::numMols] for i in range(numMols)] 

    
    ### Get relative energies associated by reference minima.
    minimaE = []
    for i, molIndices in enumerate(allIndices):
        molE = []
        for j, fileIndices in enumerate(molIndices):
            fileE = []
            for k, confNum in enumerate(fileIndices):
                if confNum == None:  #
                    fileE.append(nan)
                elif confNum == -1:   # -1 signifies reference theory
                    fileE.append(float(elists[i][j][k]))
                else:
                    fileE.append(float(elists[i][j][confNum]))
            molE.append(fileE)
        minimaE.append(molE)
    
    ### Find index for which all minima have energies (not nan) for relative Es.
    flatten = [item for sublist in minimaE for item in sublist]
    zero = 0   # initial guess
    zpass = False
    while not zpass:
        zeroths = [item[zero] for item in flatten]
        if nan in zeroths:
            zero = zero + 1
            continue
        zpass = True

    ### Take relative energies, and convert Hartrees --> kcal/mol.
    mintemp = []
    for molE in minimaE:
        temp = []
        for fileE in molE:
            temp.append([627.5095*(fileE[i]-fileE[zero]) for i in range(len(fileE))])
        mintemp.append(temp)
    minimaE = mintemp
    
    for name, rnc, minE in zip(molNames, refNumConfs, minimaE):
        PlotMolMinima(name, rnc, minE, sdfList, sdfRef, wdir, Stag=True )
        #PlotMolMinima(name, rnc, minE, sdfList, sdfRef, wdir)
