#!/usr/bin/env python

import os
import sys
import openeye.oechem as oechem
import numpy as np
from numpy import nan
import pickle
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
import procTags as pt           # for GetSDList

# Notes:
#  - if you get IndexError on line `timeSuc += `..., make sure you are
#    using the correct match.in file if you are also using debugging()

### ------------------- Functions -------------------



def compare2Mols(rmol, qmol):
    """
    For two identical molecules, with varying conformers,
        make an M by N comparison to match the M minima of
        rmol to the N minima of qmol. Match is declared
        for lowest RMSD between the two conformers and
        if the RMSD is below 0.5 Angstrom.

    Parameters
    ----------
    rmol:       reference OEChem molecule with all its filtered conformers
    qmol:       query OEChem molecule with all its filtered conformers

    """

    automorph = True   # take into acct symmetry related transformations
    heavyOnly = False  # do consider hydrogen atoms for automorphisms
    overlay = True     # find the lowest possible RMSD


    molIndices = []  # 1D list for storing indices of matched qmol confs wrt rmol

    for Rconf in rmol.GetConfs():
        print(">>>> Matching %s conformers to minima: %d <<<<"\
            % (qmol.GetTitle(),Rconf.GetIdx()+1))

        # for this Rconf, calculate/store RMSDs with all of qmol's conformers
        rsublist = []
        for Qconf in qmol.GetConfs():
            rms = oechem.OERMSD(Rconf,Qconf,automorph,heavyOnly,overlay)
            print rms
            rsublist.append(rms)

        # for this Rconf, get qmol conformer index for minimum RMSD
        thisMin=[i for i, j in enumerate(rsublist) if j == min(rsublist)][0]
        if rsublist[thisMin] <= 0.5:
            molIndices.append(thisMin)
        else:
            print('no match bc rmsd is ',rsublist[thisMin])
            molIndices.append(None)

    return molIndices


def plotMolMinima(molName, minimaE, xticklabels, selected=None,stag=False):
    # stagger plots to better see line overlap. works best with few (<4?) numFiles.
    # minimaE for ONE molecule

    refNumConfs = len(minimaE[0])
    refFile = xticklabels[0]
    numFiles = len(minimaE)

    ### Flatten this 2D list into a 1D to find min and max for plot
    flatten = [item for sublist in minimaE for item in sublist]
    floor = min(flatten)
    ceiling = max(flatten)
    if (ceiling - floor) > 4.0:
        ystep = (ceiling - floor)/9 # have 10 increments of y-axis
        ystep = round(ystep * 2) / 2 # round the step to nearest 0.5
    else:
        ystep = (ceiling - floor)

    ### Stagger each of the component files of minimaE for ease of viewing.
    if stag==True:
        tempMinimaE = []
        for i, fileE in enumerate(minimaE):
            tempMinimaE.append([x+i/2. for x in fileE])
        minimaE = tempMinimaE
        ceiling = ceiling + numFiles

    ### Figure labels.
    plttitle="Relative Energies of %s Minima" % (molName)
    plttitle+="\nby Reference File %s" % (refFile)
    ylabel="Relative energy (kcal/mol)"
    figname = "minimaE_%s.png" % (molName)

    ### Set x-axis values and labels. Can either be letter or number.

    # LETTER LABELS
    letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ' # label x-axis by letter
    rpt = (len(minimaE[0])/26)+1
    xlabs =[''.join(i) for i in itertools.product(letters, repeat=rpt)][:refNumConfs]
    # OR, NUMBER LABELS
    #xlabs = range(len(minimaE[0]))

    ### Set this up for grid.
    fig = plt.figure(figsize=(20,10))
    ax = fig.gca()
    ax.set_xticks(np.arange(-1,refNumConfs+1,2))

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle,fontsize=16)
    plt.ylabel(ylabel,fontsize=14)
    plt.xlabel("conformer minimum",fontsize=14)
    plt.xticks(range(refNumConfs),xlabs,fontsize=12)
    plt.yticks(fontsize=12)

    ### Plot the data.
    colors = mpl.cm.rainbow(np.linspace(0, 1, numFiles))
    markers = ["x","^","8","d","o","s","*","p","v","<","D","+",">","."]*10
    #for i in reversed(range(numFiles)):
    for i, FileE in enumerate(minimaE):
        if selected is not None and i not in selected:
            continue
        xi = range(refNumConfs)
        #yi = [item[i] for item in minimaE]
        yi = FileE
        plt.plot(xi,yi,color=colors[i],label=xticklabels[i],\
            marker=markers[i],markersize=9)

    ### Add legend and set plot limits.
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.xlim(-1,refNumConfs+1)
    # y axis limits: min, max, step
    ax.set_yticks(np.arange(floor-2,ceiling+2,ystep))
    plt.grid()

    plt.savefig(figname,bbox_inches='tight')
#    plt.show()
    plt.clf()

def plotAvgTimes(molName, avgTimes, sdTimes, xticklabels):
    plttitle="Conformer-Averaged Wall Times\nfor %s" % (molName)
    plttitle+="\nGeometry Optimization in Psi4"
    ylabel="time (s)"
    figname = "timebars_%s.png" % molName
    x = range(len(avgTimes))

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle,fontsize=20)
    plt.ylabel(ylabel,fontsize=18)
    plt.xticks(x,xticklabels,fontsize=14,rotation=-30, ha='left')
    plt.yticks(fontsize=14)

    ### Plot the data.
    colors = mpl.cm.rainbow(np.linspace(0, 1, len(x)))
    plt.bar(x, avgTimes, color=colors,align='center',yerr=sdTimes,ecolor='k')
    plt.savefig(figname,bbox_inches='tight')
#    plt.show()
    plt.clf()

def plotHeatRMSE(molName, rmsArray, ticklabels,ptitle='RMS error (kcal/mol)',fprefix='rmse'):
    plttitle="%s\n%s" % (ptitle,molName)
    figname = "%s_%s.png" % (fprefix, molName)
    x = range(len(rmsArray))
    y = range(len(rmsArray))

    plt.figure(figsize=(20,10))

    ### Tranpose and plot data - imshow swaps x and y
    plt.imshow(np.asarray(rmsArray).T, cmap='jet', origin='lower')
    plt.colorbar()

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle,fontsize=20)
    plt.xticks(x,ticklabels,fontsize=12,rotation=-20, ha='left')
    plt.yticks(y,ticklabels,fontsize=12)
    plt.xlabel("reference",fontsize=14)
    plt.ylabel("compared",fontsize=14)

    ### Save/show plot.
    plt.savefig(figname,bbox_inches='tight')
#    plt.show()
    plt.clf()

def plotET(molName, eneArray, timeArray, ticklabels,fprefix='scatter'):
    plttitle="RMS error vs. ratio of wall time\n%s" % molName
    figname = "%s_%s.png" % (fprefix, molName)
    colors = mpl.cm.rainbow(np.linspace(0, 1, len(eneArray)))
    markers = ["x","^","8","d","o","s","*","p","v","<","D","+",">","."]

    # use plt.plot instead of scatter to label each point
    for i, (x,y) in enumerate(zip(eneArray,timeArray)):
        plt.scatter(x,y,c=colors[i],marker=markers[i],label=ticklabels[i])

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle,fontsize=20)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    #plt.xticks(x,ticklabels,fontsize=12,rotation=-20, ha='left')
    #plt.yticks(y,ticklabels,fontsize=12)
    plt.xlabel("RMS error (kcal/mol)",fontsize=14)
    plt.ylabel("ratio of wall time",fontsize=14)

    ### Edit legend colors. All is one color since each sublist
    # colored by spectrum.
    ax = plt.gca()
    leg = ax.get_legend()
    for i in range(len(eneArray)):
        leg.legendHandles[i].set_color(colors[i])

    ### Save/show plot.
    plt.savefig(figname,bbox_inches='tight')
#    plt.show()
    plt.clf()

def shiftArray(rmsArray):
    """
    Place the first element of array in diagonal spot.
    Order of everything before and everything after is retained.
    Example,
      0 4 6                0 4 6
      0 2 7  --becomes-->  2 0 7
      0 1 3                1 3 0

    """
    for i, sublist in enumerate(rmsArray):
            sublist.insert(i,sublist.pop(0))
    return rmsArray

def matchMinima(sdfList, thryList):
    """
    For list of SDF files, match the conformer minima to those of the reference
       SDF file. Ex. Conf G of reference file matches with conf R of file3.

    Parameters
    ----------
    sdfList: str list - list of the SDF file names to be analyzed.
          This list should include reference SDF file (sdfRef) as first element.
    thryList: str list - list of levels of theory corresponding to the files
          in sdfList. E.g., ['MP2/def2-TZVP','B3LYP-D3MBJ/6-311++G**']

    Returns
    -------
    molNames: list of molecule names from the reference molecules
    refNumConfs: list of ints representing each reference mol's number of conformers
    allIndices: 2D list representing, for each sdfQuery, the conformer indices
       that match reference conformer.
       [[-1, -1, -1], [-1], [2, 3, 1], [0]] means there are two molecules, one
       with 3 confs, and one with 1 conf. The first two sublists are -1 because
       sdfQuery matched sdfRef. For the second two sublists, sdfQuery's
       mol1 conf2 matches with sdfRef's mol1 conf1.
    elists: 2D list of floats in similar format of allIndices, being energies
       of [[file1 mol1], ..., [file1 molN], [file2 mol1], ... [file2 molN]]
       Note that the mols belonging to one file are not separated in a sublist.

    """
    def loadFile(fname):
        ifs = oechem.oemolistream()
        ifs.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifs.open(fname):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % fname)
        mols = ifs.GetOEMols()
        return mols

    sdfRef = sdfList[0]
    numFiles = len(sdfList)
    allIndices = [] # for M mols, N reference minima of each mol, P matching indices for each ref minimia
    elists = [] # 2D list: K mols per file x J numFiles
    tlists = [] # 2D list: K mols per file x J numFiles
    refNumConfs = [] # number of conformers for each mol in reference file
    molNames = [] # name of each molecule. for plotting.

    for i, sdfQuery in enumerate(sdfList):
        qthry = thryList[i]
        qmethod = qthry.split('/')[0].strip()
        qbasis = qthry.split('/')[1].strip()

        print("\n\nOpening reference file %s" % sdfRef)
        molsRef = loadFile(sdfRef)

        print("Opening query file %s, and using [ %s ] energies" % (sdfQuery, qthry))
        molsQuery = loadFile(sdfQuery)

        # loop over each molecule in reference file and in query file
        for rmol in molsRef:
            #qmol = molsQuery.next()
            noMatch = True
            for qmol in molsQuery:
                if rmol.GetTitle() == qmol.GetTitle():
                    noMatch = False
                    break
            if noMatch:
                allIndices.append([-2]*rmol.NumConfs())
                elists.append([nan]*rmol.NumConfs())
                tlists.append([nan]*rmol.NumConfs())
                print("No %s molecule found in %s" % (rmol.GetTitle(), sdfQuery))
                # gotta reset the molsQuery generator
                molsQuery = loadFile(sdfQuery)
                continue

            # get energies for plotting relative energies
            elists.append(map(float, pt.GetSDList(qmol, "QM opt energy",'Psi4', qmethod, qbasis))) # adapt for SPE? === *
            tlists.append(map(float, pt.GetSDList(qmol, "opt runtime",'Psi4', qmethod, qbasis)))


            # Skip minmatch if this query file is same as reference file;
            #    before skip, get data for elists, refNumConfs, allIndices.
            if(sdfQuery == sdfRef):
                print("\nSkipping comparison against self.")
                molNames.append(rmol.GetTitle())
                refNumConfs.append(rmol.NumConfs())
                allIndices.append([-1]*rmol.NumConfs())
                continue

            # get indices of qmol conformers that match rmol conformers
            molIndices = compare2Mols(rmol, qmol)
            allIndices.append(molIndices)

    numMols = len(refNumConfs)
    molNames = molNames[:numMols]
    print "\nmolNames\n",molNames
    print "\nrefNumConfs\n",refNumConfs
    print "\nallIndices\n",allIndices
    print "\nelists\n",elists

    return molNames, refNumConfs, allIndices, elists, tlists

def getAllTimes(sdfList, thryList):
    """
    Get times saved in SD tages from files listed in python input file lines. 
       No longer used after edits to matchMinima (07/1/2017).

    Parameters
    ----------
    sdfList: str list - list of the SDF file names to be analyzed.
          This list should include reference SDF file (sdfRef) as first element.
    thryList: str list - list of levels of theory corresponding to the files
          in sdfList. E.g., ['MP2/def2-TZVP','B3LYP-D3MBJ/6-311++G**']

    Returns
    -------
    timelists: 3D list where timelists[i][j][k] is the wall time for optimizing
       for the ith level of theory, jth molecule, kth conformer

    """

    sdfRef = sdfList[0]
    timelists = []
    for i, sdfQuery in enumerate(sdfList):
        qthry = thryList[i]
        qmethod = qthry.split('/')[0].strip()
        qbasis = qthry.split('/')[1].strip()

        print("\n\nOpening reference file %s" % sdfRef)
        ifsRef = oechem.oemolistream()
        ifsRef.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifsRef.open(sdfRef):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfRef)
        molsRef = ifsRef.GetOEMols()

        print("Opening query file %s, and using [ %s ] wall times" % (sdfQuery, qthry))
        ifsQuery = oechem.oemolistream()
        ifsQuery.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifsQuery.open(sdfQuery):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfQuery)
        molsQuery = ifsQuery.GetOEMols()

        for rmol in molsRef:
            try:
                qmol = molsQuery.next()
                timelists.append(map(float, pt.GetSDList(qmol, "opt runtime",'Psi4', qmethod, qbasis))) # for opt, not spe
            except StopIteration:
                print("No %s molecule found in %s" % (rmol.GetTitle(), sdfQuery))
                timelists.append([nan]*rmol.NumConfs())
                continue

    return timelists

def calcRMSError(trimE, zeroes):
    """
    From relative energies with respect to some conformer from calcRelEne,
       calculate the root mean square error with respect to the relative
       conformer energies of the first (reference) file.

    Parameters
    ----------
    trimE: 3D list of energies, where trimE[i][j][k] represents the
      ith molecule, jth file, kth conformer rel energy
    zeroes: a 1D list of index of the reference conformer per each mol

    Returns
    -------
    relByFile: a 1D list of RMS errors for each file with reference to
      first input file

    """
    relByFile = []
    for i, molist in enumerate(trimE):
        molEnes = []
        for j, filelist in enumerate(molist):
            errs = np.asarray(filelist) - np.asarray(molist[0]) # subtract ref file
            sqrs = errs**2. # squared
            sqrs = np.delete(sqrs,zeroes[i]) # delete reference conformer (zero rel ene)
            sqrs = sqrs[~np.isnan(sqrs)] # delete nan values to get an rmse********
            mse = np.mean(sqrs)
            rmse = np.sqrt(mse)
            molEnes.append(rmse)
        relByFile.append(molEnes)

    return relByFile


def getRatioTimes(allMolTimes, zeroes):
    """
    From all molecule times, calculate relative time ratios for matched minima.
       If a conf has nan or is not matched, that time is not considered.
       After dividing by reference time, the matched conformer files are
       averaged for a particular file opt.

    Parameters
    ----------
    allMolTimes: 3D list of times, where allMolTimes[i][j][k] represents the
      ith molecule, jth file, kth conformer time (sec)

    Returns
    -------
    relByFile: a 1D list of times ratios for each file with reference to
      first input file
    sdByFile: a 1D list of standard deviation of conformer-averaged times
      relative to first input file

    """

    relByFile = []
    sdByFile = []
    for i, molist in enumerate(allMolTimes):
        molTimes = []
        molStds = []
        for j, filelist in enumerate(molist):
            rels = np.asarray(filelist)/np.asarray(molist[0])
            rels = rels[~np.isnan(rels)] # delete nan values to get avg********
            avg = np.mean(rels)
            sd = np.std(rels)
            molTimes.append(avg)
            molStds.append(sd)
        relByFile.append(molTimes)
        sdByFile.append(molStds)
    return relByFile, sdByFile

def calcRelEne(minimaE):
    """
    Calculate the relative energy. For each file, take conformer energy
       relative to minimum X. The conformer minimum is chosen from
       the first conformer for which all files have an energy value.
       Note that relative energies are taken with a file's conformers,
       not subtracting one file from another (see calcRMSError).

    Parameters
    ----------
    minimaE: 3D list of energies, where minimaE[i][j][k] represents the
      ith molecule, jth file, kth minima of that ith molecule

    Returns
    -------
    trimE: 3D list of energies as above except with relative energies
      in kcal/mol (instead of absolute energies in Hartrees). For mols
      with a single conformer it doesn't make sense to calculate rel
      energies. These mols are deleted from minimaE.
    zeroes: a 1D list of index of the reference conformer per each mol

    """

    zeroes = []
    mols2del = []
    for i, molist in enumerate(minimaE):

        # find first conformer with least nan's.
        nanCnt = []
        for j in range(len((molist[i]))):
            nanCnt.append(sum(np.isnan([item[j] for item in molist])))
        print i, nanCnt
        zeroes.append(nanCnt.index(min(nanCnt)))

#        zero = 0  # initial guess for this molecule's reference among all files
#        zpass = False  # the zero is valid if all files have it and not nan
#
#        while not zpass:
#            try: # get enes at index zero
#                print molist # infinite loop?!??!
#                zeroths = [item[zero] if sum(np.isnan(item)) < len(molist[0]) else '' for item in molist]
#            except IndexError: # error if zero index is too far
#                pass
#            if len(molist[0]) == 1: # no rel. ene for just one conf
#                mols2del.append(i)
#                print "ONLY ONE CONF FOUND FOR MOL ",i
#                break
##            elif zero >= len(molist[0]): # all values are nan
##                mols2del.append(i)
##                print "ALL NANS FOR 1+ FILE OF MOL ",i
##                zero = 0
##                break
#            elif nan in zeroths: # missing ref for 1+ file(s)
#                zero = zero + 1
#                continue
#            else: # found successful conf in all files
#                zeroes.append(zero)
#                zpass = True
#
#    # delete cases with just one conformer
#    print("ATTN: these (zero-indexed) mols were removed from analysis due to "
#         +"single conformer or no conformer matches in at least one file: ",
#         mols2del)
#    trimE = np.delete(np.asarray(minimaE),mols2del,axis=0)
#    if len(trimE) != len(zeroes):
#        print len(trimE), zeroes
#        sys.exit("Error in determining reference confs for molecules.")

    # calc relative energies, and convert Hartrees to kcal/mol.
    mintemp = []  # not sure why this is needed but writeRelEne kicks fuss without it
    for z, molE in zip(zeroes, minimaE):
    #for z, molE in zip(zeroes, trimE):
        temp = [] # temp list for this mol's relative energies
        for fileE in molE:
            temp.append([627.5095*(fileE[i]-fileE[z]) for i in range(len(fileE))])
        mintemp.append(temp)
    trimE = mintemp

    return trimE, zeroes

def writeRelEne(molName, rmse, relEnes, zero, thryList, prefix='relene'):
    """

    """

    compF = open(prefix+'_'+molName+'.dat','w')
    compF.write("# Molecule %s\n" % molName)
    compF.write("# Energies (kcal/mol) for each matched conformer relative to "+
                "conformer " + str(zero)+" across each column.\n")
    compF.write("# Rows represent conformers of this molecule; columns " +
                "represent some calculations from a particular file.\n")
    compF.write("# Columns are ordered by conformer index, then the "+
                "following levels of theory:")

    # write methods, RMSEs, integer column header
    rmsheader = "\n# "
    colheader = "\n\n# "
    for i, t in enumerate(thryList):
        compF.write("\n# %d %s" % ((i+1),t))
        rmsheader += '\t%.4f' % rmse[i]
        colheader += '\t' + str(i+1)

    compF.write("\n\n# RMS errors by level of theory, with respect to the "+
                "first level of theory listed:")
    compF.write(rmsheader)
    compF.write(colheader)

    # write each opt's relative energies
    for i in range(len(relEnes[0])):
        compF.write('\n'+str(i)+'\t')
        thisline = [x[i] for x in relEnes]
        thisline = [ '%.4f' % elem for elem in thisline ]
        thisline = '\t'.join(map(str,thisline))
        compF.write(thisline)
    compF.close()

def reorganizeSublists(theArray,allMolIndices):
    """
    Instead of grouping by files then by molecule,
        reorder to group by molecules then by file.
    Something like this:
       [[[file1 mol1] [file1 mol2]] ... [[file2 mol1] [file2 mol2]]]
    Goes to this:
       [[[file1 mol1] [file2 mol1]] ... [[file1 molN] [file2 molN]]]

    Also does checking on if minima is matched. If not, the
       value in theArray is NOT used, and nan is used instead.

    """
    minimaE = []
    for i, molIndices in enumerate(allMolIndices):
        molE = [] # all conf energies from ith mol in all files
        for j, fileIndices in enumerate(molIndices):
            fileE = []  # all conf energies from ith mol in jth file
            for k, confNum in enumerate(fileIndices):
                # None: means no conf in qmol within 0.5 Angs of rmol's conf
                # -2: means that the conformer doesn't exist (if filtered
                #   out from job not finishing, etc.)
                if confNum == None or confNum==-2:
                    fileE.append(nan)
                elif confNum == -1:   # -1 signifies reference theory
                    fileE.append(float(theArray[i][j][k]))
                else:
                    fileE.append(float(theArray[i][j][confNum]))
            molE.append(fileE)
        minimaE.append(molE)
    return minimaE

def debugging():
    molNames = ['AlkEthOH_c1178','GBI']
    refNumConfs = [19, 1]


    allIndices = [[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1], [-1], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, None, None, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 13, 14, 16, 17], [0], [0, 1, None, None, 4, 5, 6, 7, 8, None, None, 11, None, None, 14, 15, 17, 18], [0], [0, 1, None, None, 4, 5, 6, 7, 8, None, None, 11, None, None, 14, 15, 17, 18], [-2], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [-2], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0], [0, 1, 2, None, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, 14, 15, 17, 18], [0]]

    elists = [[-466.368905867, -466.36166286, -466.36349184, -466.36890724, -466.361663336, -466.361343924, -466.361344179, -466.365459408, -466.365235601, -466.363492507, -466.363195706, -466.365234831, -466.363197808, -466.361712013, -466.362098679, -466.361233835, -466.361712871, -466.362098677], [-584.717525166], [-466.180009722, -466.172137847, -466.174850662, -466.179946306, -466.172137681, -466.171684535, -466.171684698, -466.176202258, -466.175902684, -466.174850915, -466.174531697, -466.175902536, -466.174531213, -466.170415105, -466.173161792, -466.172005389, -466.179946932, -466.172677388, -466.173161702], [-584.475472907], [-465.775852253, -465.767299276, -465.770397597, -465.775652321, -465.767299753, -465.766984408, -465.76698532, -465.771940824, -465.771624392, -465.770397747, -465.769994256, -465.771628939, -465.769995123, -465.765614329, -465.768324775, -465.767058876, -465.775651356, -465.767838591, -465.76832289], [-584.013764935], [-466.195879067, -466.187978295, -466.190236868, -466.195855903, -466.187978247, -466.187365132, -466.187364931, -466.191723715, -466.191352006, -466.190237875, -466.189778842, -466.191351227, -466.189778879, -466.18598519, -466.188525019, -466.187530284, -466.195856345, -466.187942835, -466.188525035], [-584.502636745], [-465.791707241, -465.783118897, -465.785756848, -465.791562828, -465.783118227, -465.782632116, -465.782631568, -465.787431513, -465.787070752, -465.785756326, -465.785206558, -465.787071089, -465.785206994, -465.781148061, -465.783673313, -465.782575558, -465.791562833, -465.783090902, -465.78367229], [-584.040976371], [-466.204086643, -466.196172064, -466.198447456, -466.204064978, -466.196170257, -466.195575449, -466.195576595, -466.199912062, -466.199546661, -466.198447962, -466.1980118, -466.199542838, -466.198011861, -466.194150847, -466.196723667, -466.195727638, -466.204064749, -466.196159661, -466.196723726], [-584.512754158], [-466.377174577, -466.369915603, -466.371754075, -466.377148457, -466.369915597, -466.369616637, -466.36961476, -466.373708663, -466.373493021, -466.371754017, -466.371484874, -466.373491247, -466.371483481, -466.370352433, -466.369498083, -466.377147549, -466.369985869, -466.370352421], [-584.727674216], [-465.841143924, -465.833399338, -465.836198659, -465.841245956, -465.833399412, -465.833444201, -465.833444304, -465.83713488, -465.837184081, -465.836198279, -465.836263997, -465.837184139, -465.836264175, -465.832028455, -465.834381356, -465.833675455, -465.841246069, -465.834371886, -465.83438063], [-584.04531727], [-465.448164343, -465.439695447, -465.442836375, -465.448139062, -465.439695787, -465.439912519, -465.439913231, -465.444035941, -465.444126024, -465.442839778, -465.442864286, -465.444126379, -465.442863522, -465.438450484, -465.440619704, -465.439862978, -465.44814093, -465.440667701, -465.44061973], [nan], [-465.908210615, -465.90096434, -465.902602227, -465.908268604, -465.900964432, -465.900285379, -465.900284043, -465.90434806, -465.904020438, -465.902604392, -465.901900762, -465.904019637, -465.901900902, -465.898634551, -465.900967863, -465.899927061, -465.908269692, -465.900331575, -465.900967255], [-584.253752662], [-464.585265841, -464.576858317, -464.579484516, -464.585262897, -464.576857427, -464.575858466, -464.575857754, -464.580997246, -464.580466517, -464.579485523, -464.578834058, -464.580464838, -464.578834094, -464.574577198, -464.577505552, -464.575599878, -464.585262926, -464.576744295, -464.577505546], [nan], [-464.143600837, -464.134933843, -464.137334008, -464.143307447, -464.134933393, -464.134208304, -464.134207723, -464.139488586, -464.138893717, -464.137333936, -464.136607973, -464.138891561, -464.136607923, -464.131612428, -464.13509096, -464.133472079, -464.143307973, -464.134320256, -464.135091096], [-582.156839405], [-465.156073144, -465.148867482, -465.150473978, -465.155928616, -465.14886819, -465.148286217, -465.148286913, -465.15257771, -465.152203301, -465.150474003, -465.150122914, -465.152203917, -465.150122881, -465.145824728, -465.148938464, -465.147650798, -465.155928209, -465.148451922, -465.148938085], [-583.288493862], [-465.242466445, -465.235484775, -465.237129307, -465.242372664, -465.235482458, -465.2348854, -465.234885471, -465.239092284, -465.238772779, -465.237129261, -465.236817934, -465.238771756, -465.236816927, -465.232787961, -465.235675149, -465.234471301, -465.242372793, -465.235221173, -465.235672725], [-583.35286959]]
    return molNames, refNumConfs, allIndices, elists


### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
        help="Required argument on name of text file with information on\
              file(s) and levels of theory to process.\
              See README file or examples for more details. TODO")

    parser.add_argument("--readpickle", action="store_true", default=False,
        help="If specified, read in data from pickle files from each \
              directory. Input file can be same as for heat plot inputs, \
              and pickle files will be read from same directory as \
              specified output files.")

    parser.add_argument("--verbose", action="store_true", default=False,
        help="If specified, write out relative energies in kcal/mol for \
              all conformers of all mols for all files. If in doubt, \
              do specify this option.")

    parser.add_argument("--eplot",action="store_true", default=False,
        help="Generate line plots for every molecule with relative energies.")

    parser.add_argument("--tplot", action="store_true", default=False,
        help="Generate bar plots of conformer-averaged time per each \
              optimization. One plot generated per molecule.")

    parser.add_argument("--eheatplot", default=None,
        help="Specify molecule title and generate heat map of RMS errors. \
              Input file should be analogous to input file for minima matching\
               except the files should point to the relative energies .dat file\
               generated as result of minimaMatching. RMS error line is pulled\
              from each for plotting.")

    parser.add_argument("--theatplot", default=None,
        help="Specify molecule title and generate heat map of relative opt times.\
              Input file should be analogous to input file for minima matching\
               except the files should point to the relative energies .dat file\
               generated as result of minimaMatching. Data with relative times \
               are pulled from each file for plotting.")

    parser.add_argument("--etscatter", default=None,
        help="Specify molecule title and generate heat map of relative opt times.\
              Input file should be analogous to input file for minima matching\
               except the files should point to the relative energies .dat file\
               generated as result of minimaMatching.")

    args = parser.parse_args()
    opt = vars(args)
    if not os.path.exists(opt['input']):
        raise parser.error("Input file %s does not exist." % opt['filename'])
    sys.stdout.flush()


    # Read input file and store each file's information in two lists.
    sdfList = []
    thryList = []
    with open(opt['input']) as f:
        for line in f:
            if line.startswith('#'):
                continue
            dataline = [x.strip() for x in line.split(',')]
            if dataline == ['']:
                continue
            thryList.append(dataline[0])
            sdfList.append(dataline[1])


    # =========================================================================
    if not opt['readpickle']:
        molNames, refNumConfs, allIndices, elists, tlists = matchMinima(sdfList, thryList)
        pickle.dump([molNames, refNumConfs,allIndices,elists,tlists], open('match.pickle', 'wb'))
    else:
        molNames, refNumConfs, allIndices, elists, tlists = pickle.load(open('match.pickle', 'rb'))
        #molNames, refNumConfs, allIndices, elists = debugging()
    # =========================================================================

    # Reorder indices and energies lists by molecules instead of by files.
    #   [[[file1 mol1] [file2 mol1]] ... [[file1 molN] [file2 molN]]]
    #   also now the molecules are separated by sublist
    # could be done in matchMinima function but need elists for plotting
    numMols = len(refNumConfs)
    allMolIndices = [allIndices[i::numMols] for i in range(numMols)]
    elists = [elists[i::numMols] for i in range(numMols)]
    minimaE = reorganizeSublists(elists, allMolIndices)

    # =========================================================================
    trimE, zeroes = calcRelEne(minimaE)
    rmselist = calcRMSError(trimE, zeroes)

    # =========================================================================

    if opt['verbose']:
        for i, mn in enumerate(molNames):
            try:
                writeRelEne(mn, rmselist[i],trimE[i],zeroes[i],thryList)
            except IndexError:
                zeroes.append(nan)
                writeRelEne(mn, [nan]*len(thryList),elists[i],zeroes[i],thryList)


    if opt['tplot']:
#        allMolTimes = getAllTimes(sdfList, thryList) # ordered by file, mol, conf
        allMolTimes = tlists

        # match conformer times using indices from matchMinima then get stats
        allFileTimes = [[] for i in range(numMols)]
        allFileStds = [[] for i in range(numMols)]
        for i in range(len(sdfList)*numMols):
            timeSuc = 0     # running sum of successfully matched minima times
            numSuc = 0      # running number of successful minima
            numconfs = refNumConfs[i%numMols]  # i%numMols gets some mol
            fileTimes = []  # collect successful times for stdevs
            for k in range(numconfs):
                thisIndex = allIndices[i][k]
                if thisIndex == -1:
                    timeSuc += allMolTimes[i][k]
                    numSuc += 1
                    fileTimes.append(allMolTimes[i][k])
                elif thisIndex > -1:
                    timeSuc += allMolTimes[i][thisIndex]
                    numSuc += 1
                    fileTimes.append(allMolTimes[i][thisIndex])
            try:
                fTimeAvg = float(timeSuc)/numSuc
            except ZeroDivisionError:
                fTimeAvg = nan
            allFileTimes[i%numMols].append(fTimeAvg)
            allFileStds[i%numMols].append(np.std(np.array(fileTimes)))
            #print 'average ',fTimeAvg,' over ',numSuc," samples"

        # separately, go from allMolTimes and calculate relative speeds
        allMolTimes = [allMolTimes[i::numMols] for i in range(numMols)]
        timesByMol = reorganizeSublists(allMolTimes, allMolIndices)
        relTimes, sdTimes = getRatioTimes(timesByMol, zeroes)

        # bar plot of average times with stdevs
        for name, fileTimes, stdevs in zip(molNames, allFileTimes, allFileStds):
            print name, 'times ',fileTimes
            print name, 'stdevs ',stdevs
            plotAvgTimes(name, fileTimes, stdevs, thryList)
#            plotAvgTimes(name, fileTimes[1:], stdevs[1:], thryList[1:]) # ================== * REMOVE last element from all

        if opt['verbose']: # append time to relative energies file
            for i, name in enumerate(molNames):
#            for name, fileTimes, stdevs in zip(molNames, allFileTimes, allFileStds, relTimes, sdTimes):
                compF = open('relene_'+name+'.dat','a')
                compF.write("\n\n# Avg times, stdevs, avg time ratios relative to ref, stdev of rel. time ratios:")
                avgline = "\n# "
                stdline = "\n# "
                a2line =  "\n# "
                s2line =  "\n# "
                for j, t in enumerate(thryList):
                    avgline += ' %.4f' % allFileTimes[i][j]
                    stdline += ' %.4f' % allFileStds[i][j]
                    a2line += ' %.4f' % relTimes[i][j]
                    s2line += ' %.4f' % sdTimes[i][j]

                compF.write(avgline)
                compF.write(stdline)
                compF.write(a2line)
                compF.write(s2line)
                compF.close()


    if opt['eplot']:
        for name, minE in zip(molNames, trimE):
            plotMolMinima(name, minE, thryList)
            #plotMolMinima(name, minE, thryList, selected=[0,7,12]) # zero based index

    if opt['eheatplot'] is not None:
        rmsArray = []
        for infile in sdfList:
            with open(infile) as f:
                for line in f:
                    if "RMS error" in line:
                        rmse = itertools.islice(f,1).next()
                        rmse = [float(s) for s in rmse.split()[1:]]
                        rmsArray.append(rmse)
                        break
        rmsArray = shiftArray(rmsArray)
        for t in rmsArray:
            print t
        plotHeatRMSE(opt['eheatplot'],rmsArray,thryList)

    if opt['theatplot'] is not None:
        rmsArray = []
        for infile in sdfList:
            with open(infile) as f:
                for line in f:
                    if "avg time" in line:
                        ravgs = itertools.islice(f,3) # get line via iterator
                        for j in ravgs:            # get last item of iterator
                            pass
                        ravgs = [float(s) for s in j.split()[1:]]
                        rmsArray.append(ravgs)
                        break
        rmsArray = shiftArray(rmsArray)
#        rmsArray = [item[:-1] for item in rmsArray] # ================== * REMOVE last element from all
        plotHeatRMSE(opt['theatplot'],rmsArray,thryList, ptitle='Ratio of wall times',fprefix='times')

    if opt['etscatter'] is not None:
        eArray = []
        tArray = []
        for infile in sdfList:
            with open(infile) as f:
                for line in f:
                    if "RMS error" in line:
                        rmse = itertools.islice(f,1).next()
                        rmse = [float(s) for s in rmse.split()[1:]]
                        eArray.append(rmse)
                    if "avg time" in line:
                        ravgs = itertools.islice(f,3) # get line via iterator
                        for j in ravgs:            # get last item of iterator
                            pass
                        ravgs = [float(s) for s in j.split()[1:]]
                        tArray.append(ravgs)
                        break
        eArray = shiftArray(eArray)
        tArray = shiftArray(tArray)
        eArray = [item[:-1] for item in eArray] # ================== * REMOVE last element from all
        tArray = [item[:-1] for item in tArray] # ================== * REMOVE last element from all
        print len(eArray)
        for i in range(len(sdfList)):
#            if i<13: continue
            plotET(opt['etscatter'],eArray[i],tArray[i],thryList,fprefix='scatter'+str(i+1))
