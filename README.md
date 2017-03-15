
# Overview

# Components

  * 12\_executor.py
  * confs2psi.py
  * examples
  * filterConfs.py
  * getPsiResults.py
  * matchMinima.py
  * selectConfs.tcl
  * smi2confs.py
  * timeAvg.py
  * viewer.ipynb
 
# Naming

Base names (e.g. basename.smi, basename.sdf) can contain underscores but no dashes or dots.
  * Dash is used for SDF numbering code, dot is used for splitting based on file extension.
  * Use an underscore if you want.

Smiles file should contain, in each line: "SMILESSTRING title" and be named in format of filename.smi.
  * Molecule title should have no dashes, as Psi4 will raise an error.

Examples:
  * CC(C(C(C)O)O)O AlkEthOH\_c42
  * CCCC AlkEthOH\_c1008
  * CCOC(C)(C)C(C)(C)O AlkEthOH\_c1178

# Instructions .... update me!
Execute these commands in the directory that you want input/output files to be generated.

 1. Generate conformers, quick MM optimization, Psi4 input files.
     * python executor.py -f /include/full/path/to/file.smi --setup -m 'mp2' -b 'def2-sv(p)'

 2. Get Psi4 results from the last set of optimizations.
     * python executor.py -f /include/full/path/to/file-200.sdf --results -m 'mp2' -b 'def2-sv(p)'

 3. In a new (sub?)directory, set up Psi4 SPE calculations from last results.
     * python executor.py -f /include/full/path/to/file-220.sdf --setup --spe -m 'b3lyp-d3mbj' -b 'def2-tzvp'

 4. Get Psi4 results from SPEs.
     * python executor.py -f /include/full/path/to/file-220.sdf --results --spe -m 'b3lyp-d3mbj' -b 'def2-tzvp'

 5. Join results of SPEs and plot energy values (not recd for test sets > 20). Else can plot on a single mol's conformers.
     * python executor.py -f       --combine

 6. Compare SPE and OPT2 files.
     * Copy both files in new directory (optional).
     * Open python and "import timeAvg.py". Must be in same directory as script.
     * timeAvg.compareSPEopt('/include/path/to/file1.ext','/include/path/to/file2.ext', "QM spe", "QM opt energy", 'b3lyp-d3mbj','def2-tzvp',m2=None,b2=None,verbose=True)

# Output

 SDF file numbering codes (along the lines of UNIX permissions...)
 Here, x is used as a placeholder for either 0, 1, or 2.

 * 000 original file
 * 1xx MM opt but no filter
 * 2xx MM opt and filter
 * x1x QM opt but no filter
 * x2x QM opt and filter
 * xx1 either QM second opt or SPE and no filter
 * xx2 either QM second opt or SPE and filter

f prefix means filtered but no MM. So like f020 means filtered from orig, no MM opt/filter, QM opt/filter.
I'm not making this its own category since MM opt *should* be in the regular workflow.

So the starting code is basename-000.sdf and the max code is basename-222.sdf.
Particularly relevant files are 000, 200, 220, 222. 

