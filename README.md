
# Overview

# Components

  * avgTimeEne.py
  * confs2psi.py
  * examples/
  * executor.py
  * filterConfs.py
  * getPsiResults.py
  * selectConfs.tcl
  * smi2confs.py
  * stitchSpe.py
  * viewer.ipynb
  * writeOneConf.py
 

# Instructions .... update me!
Execute these commands in the directory that you want input/output files to be generated.
Before starting, you need an input file with a list of SMILES strings and corresponding molecule titles.
See section on "Naming molecules in the input SMILES file" and "File name limitations".

 1. Generate conformers, quick MM optimization, Psi4 input files.
     * python executor.py -f /include/full/path/to/file.smi --setup -m 'mp2' -b 'def2-sv(p)'

 2. Get Psi4 results from the last set of optimizations.
     * python executor.py -f /include/full/path/to/file-200.sdf --results -m 'mp2' -b 'def2-sv(p)'

 3. In a new (sub?)directory, set up Psi4 SPE calculations from last results.
     * python executor.py -f /include/full/path/to/file-220.sdf --setup --spe -m 'b3lyp-d3mbj' -b 'def2-tzvp'

 4. Get Psi4 results from SPE.
     * python executor.py -f /include/full/path/to/file-220.sdf --results --spe -m 'b3lyp-d3mbj' -b 'def2-tzvp'

 5. Combine results from various job types to calculate model uncertainty.
     * See section on "Creating input file for stitchSpe.py"
     * mpython /data12/cmf/limvt/qm_AlkEthOH/pipeline/01_scripts/stitchSpe.py -i /path/and/input.dat --barplots

 6. (opt.) If some mol has a high RMSD, identify conformer and visualize structure.
     * mpython /data12/cmf/limvt/qm_AlkEthOH/pipeline/01_scripts/writeOneConf.py ...............

 7. (opt.) Get wall clock times, num opt steps, relative energies. 
     * mpython /data12/cmf/limvt/qm_AlkEthOH/pipeline/01_scripts/avgTimeEne.py --relene -f /path/&/file.sdf -m 'b3lyp-d3mbj' -b 'def2-tzvp'



# Naming molecules in the input SMILES file

Smiles file should contain, in each line: "SMILESSTRING title" and be named in format of filename.smi.
  * Molecule title should have no dashes, as Psi4 will raise an error.

Examples:
  * CC(C(C(C)O)O)O AlkEthOH\_c42
  * CCCC AlkEthOH\_c1008
  * CCOC(C)(C)C(C)(C)O AlkEthOH\_c1178



# File name limitations

Base names (e.g. basename.smi, basename.sdf) can contain underscores but NO dashes or dots.
  * Dash is used for SDF numbering code (see below).
  * Dot is used for splitting based on file extension.


# Output files throughout the pipeline

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
Or the 221 files for comparing relative energies of single pt energy calcns.

# Creating input file for `stitchSpe.py`

 * This should be a text file directing the script to process a particular quantity.
 * The first uncommented line should be the keyword of the specific quantity (e.g., energy) found in the SD tag label.
 * Following lines should contain the following information in order, separated by a comma:
    * sdf file with full path
    * Boolean, True if SPE, False if optimization
    * method
    * basis set
 * The first sdf file listed will be the reference values for all following lines when computing RMSDs.
 * The sdf files on each line should ALL have the same molecules, same conformers, etc. These may differ in coordinates or SD tags.

 * Example:
    # comments begin with pound symbol and are ignored
    energy
    /path/and/setofMols-221-opt2.sdf, False, b3lyp-d3mbj ,    def2-tzvp   
    /path/and/setofMols-221-spe1.sdf, True , b3lyp-d3mbj ,    def2-tzvp   
    /path/and/setofMols-221-spe2.sdf, True , mp2, cc-pvtz ,  
    /path/and/setofMols-221-spe3.sdf, True , pbe0, 6-311g**

