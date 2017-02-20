#!/usr/bin/env python

import os, sys
import argparse

import smi2confs
import filterConfs
import confs2psi
import getPsiResults

# TODO
# make method and basis set command line arguments bc opt2 will differ

# Usage: python executor.py -f /path/to/inputfile --setup 

# Envisioned pipeline stages:
#   I.   feed in smiles file > generate confs > MM optimize > generate psi4 files
#   II. 
#   III. feed in sdf file with QM coords > generate SPE psi4 input files 

# Flexibility: 
#   - can feed it in an SDF file (which already has ready-2-go-confs) to set up Psi4.
#   - in special cases, can feed in a mol2 file for conformers of ONE molecule.
#     User will have to specify molecule name and charge (for ALL) confs.
#     If these differ among confs, don't go the mol2 route. Test case make all charges -6
#     instead of neutral. Mol name may be ok, but used VMD which reset names.


# Note 1: This pipeline uses some preset parameters, such as 
#    resClash=True and quickOpt=True (with SD opt) in smi2confs, and
#    MP2/def2-sv(p) for QM opt1. These can be modified in the argument
#    inputs here, or in the parent code itself.
# Note 2: The input file must be in the same directory that the script is
#    called. The input files are generated in a subdir in this dir.


def main(**kwargs):
    _, extension = os.path.splitext(opt['filename'])
    adir, fname = os.path.split(opt['filename'])
    if adir == '': adir = './'
    base = fname.replace('-', '.').split('.')[0]

    if opt['setup']:

        if extension == '.smi': ### MM opt & filter
            print("\nGenerating and filtering conformers for %s in %s" % (opt['filename'], adir))
            msdf = base + '.sdf'
            smi2confs.smi2confs(adir, opt['filename'])
            filterConfs.filterConfs(adir, msdf, "MM Szybki SD Energy", 200)
        else: msdf = opt['filename']

        ### Generate Psi4 inputs.
        print("\nCreating Psi4 input files for %s in %s ..." % (base, adir))
        if not opt['spe']:
#            confs2psi.confs2psi(adir,msdf,'b3lyp-d3mbj','def2-tzvp',False,"1.5 Gb")
            confs2psi.confs2psi(adir,msdf,'mp2','def2-sv(p)',False,"1.5 Gb")
        else:
            confs2psi.confs2psi(adir,msdf,'b3lyp-d3mbj','def2-tzvp',True,"1.5 Gb")

    else:  # ========== AFTER QM =========== #
        print("Getting Psi4 results and filtering for %s ..." %(opt['filename']))

        ### Specify output file name
        if "220" not in opt['filename']:
            osdf = base + '-210.sdf'
            suffix = '220'
        else:
            osdf = base + '-221.sdf'
            suffix = '222'
        if os.path.exists(osdf):
            print("File already exists: %s. Exiting without getting results.\n" % (osdf))
            return

        ### Get results and filter
        method, basisset = getPsiResults.getPsiResults(adir, opt['filename'], osdf, spe=opt['spe'])
        if not opt['spe']: tag = "QM Psi4 Final Opt. Energy (Har) b3lyp-d3mbj/def2-tzvp"
        else: tag = "QM Psi4 Single Pt. Energy (Har) %s/%s" % (method, basisset)
        filterConfs.filterConfs(adir, osdf, tag, suffix)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    req = parser.add_argument_group('required arguments')

    req.add_argument("-f", "--filename",
        help="SDF file (with FULL path) to be set up or processed.")

    parser.add_argument("--setup", action="store_true", default=False,
        help="If True (default=False), generate and filter conformers (if\
 starting from *.smi), then generate Psi4 input files.")

    parser.add_argument("--results", action="store_true", default=False,
        help="If True (default=False), process Psi4 output files\
 and filter conformers.")

    parser.add_argument("--spe", action="store_true", default=False,
        help="If True, either set up or process single point energy\
 calculations. If False, will set up or process geometry optimizations.\
 Can be used with either the --setup flag or the --results flag.")

    args = parser.parse_args()
    opt = vars(args)

    ### Check that both 'setup' and 'spe' are not both true or both false.
    if opt['setup'] == opt['results']:
        raise parser.error("Exactly one of either --setup or --results must\
 be specified.")

    ### Check that input file exists.
    if not os.path.exists(opt['filename']):
        raise parser.error("Input file %s does not exist. Try again." % opt['filename'])

    ### Check that all dependencies are present.
#    if opt['setup']:
#        if opt['smiles'] == None:
#            raise argparse.ArgumentError("No SMILES file provided for setup.")


    main(**opt)

