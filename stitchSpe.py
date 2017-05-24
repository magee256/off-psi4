
## DONE: make so that cc-pVTZ == cc-pvtz == whatever else
## DONE: edit writer so that it goes in main, and format file
##   to have file header, mol header, mol columns (conf, ref, spe1, spe2, ...)
## DONE: plot: read in mol title, plot like with other case, large dots
## TODO: check to make sure all confNumsare the same, at leats the same length?

# Note: If you see error: "ValueError: could not convert string to float:"
#   check to make sure that the value for SDF tag is correct.
import os
import openeye.oechem as oechem
import numpy as np
import argparse
import procTags as pt
import collections # ordered dictionary

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator as o

def barplot(ax, dpoints):
    '''
    This function written by Peter Kerpedjiev.
    http://emptypipes.org/2013/11/09/matplotlib-multicategory-barchart/

    Create a barchart for data across different categories with
    multiple conditions for each category.

    @param ax: The plotting axes from matplotlib.
    @param dpoints: The data set as an (n, 3) numpy array
    '''

    # Aggregate the conditions and the categories according to their
    # mean values
    conditions = [(c, np.mean(dpoints[dpoints[:,0] == c][:,2].astype(float)))
                  for c in np.unique(dpoints[:,0])]
    categories = [(c, np.mean(dpoints[dpoints[:,1] == c][:,2].astype(float)))
                  for c in np.unique(dpoints[:,1])]

    # sort the conditions, categories and data so that the bars in
    # the plot will be ordered by category and condition
    conditions = [c[0] for c in sorted(conditions, key=o.itemgetter(1))]
    categories = [c[0] for c in sorted(categories, key=o.itemgetter(1))]

    dpoints = np.array(sorted(dpoints, key=lambda x: categories.index(x[1])))

    # the space between each set of bars
    space = 0.3
    n = len(conditions)
    width = (1 - space) / (len(conditions))

    # Create a set of bars at each position
    for i,cond in enumerate(conditions):
        indeces = range(len(categories))
        #indeces = range(1, len(categories)+1)
        vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.float)
        pos = [j - (1 - space) / 2. + i * width for j in indeces]
        ax.bar(pos, vals, width=width, label=cond,
               color=cm.Accent(float(i) / n))

    # Set the x-axis tick labels to be equal to the categories
    ax.set_xticks(indeces)
    ax.set_xticklabels(categories)
    plt.setp(plt.xticks()[1], rotation=50)

    # Add the axis labels
    ax.set_ylabel("energy (kcal/mol)")
    ax.set_title("RMSDs of relative conformer energies")

    # Add a legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left')


def calcRelEne(dict1, dict2):
    """

    WORDS

    Parameters--------UPDATE ME
    ----------
    sdf1 | str | path+name of SDF file with times for all confs of all mols
    sdf2 | str | should have same mols/confs as sdf1, likely diff coords/tags
    tag1 | str | data in this tag from sdf1 to be compared to sdf2
    tag2 | str | data in this tag from sdf2 to be compared to sdf1
    m1/b1| str | method/basis from sdf1. If m2/b2 is None, use same from m1/b1.

    For tags, see options in GetSDList function.

    """

    def prelim(sdfRef,spe):
        # Open file.
        ifs1 = oechem.oemolistream()
        ifs1.SetConfTest( oechem.OEAbsoluteConfTest() )
        if not ifs1.open(sdfRef):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfRef)
        mols = ifs1.GetOEMols()

        # Determine SD tag from which to obtain energy.
        if spe.lower()=='true': tagword = "QM spe"
        else: tagword = "QM opt energy"
        return mols, tagword

    mols1, tag1 = prelim(dict1['fname'],dict1['fromspe'])
    mols2, tag2 = prelim(dict2['fname'],dict2['fromspe'])
    headerMols = []
    titleMols = []
    rmsds = []
    confNums = []
    refEnes = []
    compEnes = []

    for imol in mols1:
        jmol = mols2.next()

        # Get absolute energies from the SD tags
        #print dict2['fromspe'],tag2, dict2['method'],dict2['basisset'] # for debugging
        #print pt.GetSDList(jmol, tag2, dict2['method'],dict2['basisset']) # for debugging
        iabs = np.array(map(float, pt.GetSDList(imol, tag1, dict1['method'],dict1['basisset'])))
        jabs = np.array(map(float, pt.GetSDList(jmol, tag2, dict2['method'],dict2['basisset'])))

        # Get omega conformer number of first, for reference info
        # whole list can be used for matching purposes
        originum = pt.GetSDList(imol, "original index")
        origjnum = pt.GetSDList(jmol, "original index")


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

        # convert relative energies from Hartrees to kcal/mol
        irel = 627.5095*irel
        jrel = 627.5095*jrel

        header = ("\n# Mol %s, RMSD = %.5f kcal/mol\n" % (imol.GetTitle(), rt))
        header += ("# Energies relative to omega conf #%s\n" % (originum[0]))
        headerMols.append(header)
        titleMols.append(imol.GetTitle())
        rmsds.append(rt)
        confNums.append(originum)
        refEnes.append(irel)
        compEnes.append(jrel)

    return headerMols, titleMols, rmsds, confNums, refEnes, compEnes



### ------------------- Script -------------------

def main(wholedict, verbose=False,outfn='relene-rmsd.dat', plotbars=False):
    """
    Parameters
    ----------
    wholedict
    """

    sdfRef = wholedict[0]['fname']
    print("Using reference file: %s " % sdfRef)

    # Write description in output file.
    compF = open(os.path.join(os.path.dirname(wholedict[0]['fname']),outfn), 'w')
    compF.write("# RMSD of relative energies (kcal/mol):\n")

    for i, d in enumerate(wholedict.values()):
        compF.write("# File %d: %s\n" % (i, d['fname']))
        compF.write("#   using %s/%s energy, SPE=%s\n" % (d['method'],d['basisset'], d['fromspe']))

    for i in range(1,len(wholedict)):
        compfile = wholedict[i]['fname']
        print("Starting comparison on file: %s" % compfile)

        # each of the four returned vars is (file) list of (mols) lists
        headerMols, titleMols, rmsds,confNums, refEnes, compEnes =\
               calcRelEne(wholedict[0], wholedict[i])
        wholedict[i]['headerMols'] = headerMols
        wholedict[i]['titleMols'] = titleMols
        wholedict[i]['rmsds'] = rmsds
        wholedict[i]['confNums'] = confNums
        wholedict[i]['refEnes'] = refEnes
        wholedict[i]['compEnes'] = compEnes

    # loop over each mol and write energies from wholedict by column
    for m in range(len(wholedict[1]['titleMols'])):
        compF.write('\n\n# Mol '+wholedict[1]['titleMols'][m])

        # for this mol, write the rmsd from each file side by side
        line = ' RMSDs: '
        for i in range(1,len(wholedict)):
            line += (str(wholedict[i]['rmsds'][m])+'\t')
        compF.write(line)
        if not verbose:
            continue

        compF.write('\n# ==================================================================')
        compF.write('\n# conf\trel. enes. in column order by file (listed at top)')


        # for this mol, write the compEnes from each file by columns
        for c in range(len(wholedict[1]['confNums'][m])):
            line = '\n'+str(wholedict[1]['confNums'][m][c])+'\t'
            line += str(wholedict[1]['refEnes'][m][c])+'\t'
            for i in range(1,len(wholedict)):
                line += (str(wholedict[i]['compEnes'][m][c])+'\t')
            compF.write(line)

    compF.close()


    if plotbars:
        # fill in dictionary for ref file for list comprehension later
        wholedict[0]['titleMols'] = np.full(len(wholedict[1]['titleMols']), np.nan)
        wholedict[0]['rmsds'] = np.full(len(wholedict[1]['rmsds']), np.nan)

        # extract part of dictionary using list comprehension
        subset = np.array([[(wholedict[fd][key]) for key in\
('ftitle','titleMols','rmsds')] for fd in wholedict.keys()[1:]], dtype=object).T

        # build plot list
        plotlist = []
        for m in range(len(wholedict[1]['titleMols'])):
            for f in range(len(wholedict)-1):
                temp = []
                #temp.append(subset[0][f].split('/')[-1].split('.')[0])
                temp.append(subset[0][f])
                temp.append(subset[1][f][m])
                temp.append(subset[2][f][m])
                plotlist.append(temp)
        plotlist = np.array(plotlist)
    
        # generate plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        barplot(ax, plotlist)
        plt.savefig('barchart.png',bbox_inches='tight')
        plt.show()

    #print np.array([[(wholedict[fd]['fname'],wholedict[fd][key]) for key in ('titleMols','rmsds')] for fd in wholedict.keys()[1:]], dtype=object)
    #print np.array([[wholedict[fd][key] for key in ('titleMols','rmsds')] for fd in wholedict.keys()[1:]]).T
    #print np.array([(wholedict[fd]['rmsds']) for fd in wholedict.keys()])
    #print np.array([(wholedict[fd]['titleMols'], wholedict[fd]['rmsds']) for fd in wholedict.keys()])
    return



### ------------------- Parser -------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
        help="Required argument on name of text file with information on\
              file(s) and levels of theory to process.\
              See README file or examples for more details. TODO")

    parser.add_argument("--verbose", action="store_true", default=False,
        help="If specified, write out relative energies in kcal/mol for \
              all conformers of all mols for all files.")

    # TODO
    parser.add_argument("-l", "--lineplot",
        help="Optional argument, name of text file with molecule name(s)\
              for which line plots will be generated of relative energies\
              of all compared quantities.\
              See README file or examples for more details. TODO")

    parser.add_argument("--plotbars", action="store_true", default=False,
        help="If specified, generate bar plots of each comparison file wrt \
              reference file (first entry of input file).")

    args = parser.parse_args()
    opt = vars(args)
    if not os.path.exists(opt['input']):
        raise parser.error("Input file %s does not exist." % opt['filename'])

    # Read input file and store each file's information in an overarching set.
    # http://stackoverflow.com/questions/25924244/creating-2d-dictionary-in-python
    linecount = 0
    wholedict = collections.OrderedDict()
    with open(opt['input']) as f:
        for line in f:
            if line.startswith('#'):
                continue
            dataline = [x.strip() for x in line.split(',')]
            wholedict[linecount] = {'ftitle':dataline[0],'fname':dataline[1], 'fromspe':dataline[2], 'method':dataline[3], 'basisset':dataline[4]}
            linecount += 1

    main(wholedict, opt['verbose'],plotbars=opt['plotbars'])
