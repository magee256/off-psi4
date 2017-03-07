
# Overview

# Components
 
# Naming

 Base names (e.g. basename.smi, basename.sdf) can contain underscores but no dashes or dots.
  > Dash is used for SDF numbering code, dot is used for splitting based on file extension.
  > Use an underscore if you want.

 Smiles file should contain, in each line: "SMILESSTRING title" and be named in format of filename.smi.
  > Molecule title should have no dashes, as Psi4 will raise an error.
 Example:
    CC(C(C(C)O)O)O AlkEthOH_c42
    CCCC AlkEthOH_c1008
    CCOC(C)(C)C(C)(C)O AlkEthOH_c1178

# Instructions .... update me!

 xx. Generate conformers, quick MM optimization, Psi4 input files.
     python /work/cluster/limvt/qm_AlkEthOH/pipeline/01_scripts/12_executor.py -f /include/full/path/to/file.smi --setup

 xx. Get Psi4 results from the last set of optimizations.
     python /work/cluster/limvt/qm_AlkEthOH/pipeline/01_scripts/12_executor.py -f /include/full/path/to/file-200.sdf --results

 xx. Set up Psi4 SPE calculations from last results.
     python /work/cluster/limvt/qm_AlkEthOH/pipeline/01_scripts/12_executor.py -f /include/full/path/to/file-220.sdf --setup --spe

 xx. Get Psi4 results from SPEs.
     python /work/cluster/limvt/qm_AlkEthOH/pipeline/01_scripts/12_executor.py -f /include/full/path/to/file-220.sdf --results --spe

 xx. Compare SPE and OPT2 files.
     Copy both files in new directory (optional).
     Open python and "import timeAvg.py". Must be in same directory as script.
     timeAvg.compareSPEopt('/include/path/to/file1.ext','/include/path/to/file2.ext', "QM spe", "QM opt energy", 'b3lyp-d3mbj','def2-tzvp',m2=None,b2=None,verbose=True)

# Output

 SDF file numbering codes (along the lines of UNIX permissions...)
 Here, x is used as a placeholder for either 0, 1, or 2.

  -000 original file
  -1xx MM opt but no filter
  -2xx MM opt and filter
  -x1x QM opt but no filter
  -x2x QM opt and filter
  -xx1 either QM second opt or SPE and no filter
  -xx2 either QM second opt or SPE and filter

 f prefix means filtered but no MM. So like f020 means filtered from orig, no MM opt/filter, QM opt/filter.
 I'm not making this its own category since MM opt *should* be in the regular workflow.

 So the starting code is basename-000.sdf and the max code is basename-222.sdf.
 Particularly relevant files are 000, 200, 220, 222. 

