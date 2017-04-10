#!/usr/bin/env python

import os, sys
import argparse

import smi2confs
import filterConfs
import confs2psi
import getPsiResults

# Envisioned pipeline stages with this script:
#   I.   feed in smiles file > generate confs > MM optimize > generate psi4 inputs
#   II.  process psi4 results
#   III. generate new Psi4 inputs from II (e.g. SPE or OPT2)
#   IV.  process psi4 results

# Example usage:
#   python executor.py -f /path/to/inputfile --setup -m 'mp2' -b 'def2-sv(p)'
#   python executor.py -f /path/to/inputfile --results -m 'mp2' -b 'def2-sv(p)'
#   python executor.py -f /path/to/inputfile --setup --spe -m 'b3lyp-d3mbj' -b 'def2-tzvp'
#   python executor.py -f /path/to/inputfile --results --spe -m 'b3lyp-d3mbj' -b 'def2-tzvp'

# Note 1: This pipeline uses some preset parameters, such as 
#    resClash=True and quickOpt=True (with SD opt) in smi2confs, and
#    MP2/def2-sv(p) for QM opt1. These can be modified in the argument
#    inputs here, or in the parent code itself.
# Note 2: The input file must be in the same directory that the script is
#    called. The directory tree goes (pwd)/molName/confNum .

# Flexibility: 
#   - can feed it in an SDF file (which already has ready-2-go-confs) to set up Psi4.
#   - can sort of setup mol2 files (e.g. for one mol and all its confs) but 
#     check that molecule name and total charge is correct in Psi4 input files.

def main(**kwargs):
    _, extension = os.path.splitext(opt['filename'])
    adir, fname = os.path.split(opt['filename'])
    if adir == '' or adir is None or adir is '.':
        adir = os.getcwd()
        fullname = os.path.join(adir,opt['filename'])
    else:
        fullname = opt['filename']
    base = fname.replace('-', '.').split('.')[0]

    if opt['setup']:

        if extension == '.smi': ### MM opt & filter
            print("\nGenerating and filtering conformers for %s in %s" % (opt['filename'], adir))
            msdf = base + '.sdf'
            smi2confs.smi2confs(adir, opt['filename'])
            filterConfs.filterConfs(adir, msdf, "MM Szybki SD Energy", 200)
        else: msdf = fullname

        ### Generate Psi4 inputs.
        print("\nCreating Psi4 input files for %s in %s ..." % (base, adir))
        confs2psi.confs2psi(msdf,opt['method'],opt['basisset'],opt['spe'],"5.0 Gb")

    else:  # ========== AFTER QM =========== #

        ### Specify output file name
        if "220" not in fname:
            osdf = base + '-210.sdf'
            suffix = '220'
        else:
            osdf = base + '-221.sdf'
            suffix = '222'

        ### Get results.
        print("Getting Psi4 results for %s ..." %(fname))
        method, basisset = getPsiResults.getPsiResults(fullname, osdf, spe=opt['spe'])
        if method is None or basisset is None:
            method = opt['method']
            basisset = opt['basisset']

        ### Filter.
        print("Filtering Psi4 results for %s ..." %(fname))
        if not opt['spe']: 
            tag = "QM Psi4 Final Opt. Energy (Har) %s/%s" % (method, basisset)
        else: 
            tag = "QM Psi4 Single Pt. Energy (Har) %s/%s" % (method, basisset)
        filterConfs.filterConfs(os.path.join(adir,osdf), tag, suffix)



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

    req.add_argument("-m", "--method",
        help="Name of QM method. Put this in 'quotes'.")
    req.add_argument("-b", "--basisset",
        help="Name of QM basis set. Put this in 'quotes'.")

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

