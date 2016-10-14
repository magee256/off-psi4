#!/usr/bin/env python

import os, sys
import smi2confs
import filterConfs
import confs2psi
import getPsiResults
import matchMinima

# directories should end in /
# sdf file should have no . except for extension

adirs = []
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-A/' ]  
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-A/'  ] 
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-A/']  
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-A/']  

bdirs = []
bdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-B/' ]  
bdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-B/'  ] 
bdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-B/']  
bdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-B/']  

cdirs = []
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-C/' ]  
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-C/'  ] 
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-C/']  
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-C/']  

ddirs = []
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-D/' ]  
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-D/'  ] 
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-D/']  
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-D/']  

smiles = []
smiles += ['AlkEthOH_chain_tiny.smi' ] 
smiles += ['AlkEthOH_rings_tiny.smi' ] 
smiles += ['AlkEthOH_rings_supp1.smi'] 
smiles += ['diverse.smi'             ] 


### For A, C directories: getPsiResults, filter (after QM).
for i, smi in enumerate(smiles):
    sdf = smi.split('.')[0] + '-minima.sdf'
    adir = adirs[i]
    cdir = cdirs[i]
    osdf = smi.split('.')[0] + '-optimized.sdf'

    if i==1 or i==2: # only do comparison for C-1 and C-3
        continue

#    print("Getting Psi4 results for file %d: %s" %(i+1, smi))
#    getPsiResults.getPsiResults(adir, sdf, osdf)
#    getPsiResults.getPsiResults(cdir, sdf, osdf)
#
#    print("Filtering SDF molecules for file %d: %s" %(i+1, smi))
#    filterConfs.filterConfs(adir, osdf, "QM Psi4 Final Opt. Energy (Har) mp2/def2-sv(p)")
#    filterConfs.filterConfs(cdir, osdf, "QM Psi4 Final Opt. Energy (Har) mp2/def2-sv(p)")

for smi, adir, bdir in zip(smiles, adirs, bdirs):
     base = smi.split('.')[0]
#     confs2psi.prep(adir,adir+'SPE/', base+'-optimized-minima.sdf', base+'-SPEpre.sdf')
     confs2psi.prep(adir,bdir+'OPT2/', base+'-optimized-minima.sdf', base+'-OPT2pre.sdf')

#     confs2psi.confs2psi(adir+'SPE/', base+'-SPEpre.sdf', 'b3lyp-d3mbj','def2-tzvp', True, "1.5 Gb") 
     confs2psi.confs2psi(bdir+'OPT2/', base+'-OPT2pre.sdf', 'b3lyp-d3mbj','def2-tzvp', False, "1.5 Gb") 


#### For the D directories, confs2psi with green star method.
#for wdir, smi in zip(ddirs, smiles):
#    sdf = smi.split('.')[0] + '-minima.sdf'
#
#    ### After QM
#    osdf = smi.split('.')[0] + '-optimized.sdf'
#    getPsiResults.getPsiResults(wdir, sdf, osdf)
#    filterConfs.filterConfs(wdir, osdf, "QM Psi4 Final Opt. Energy (Har) b3lyp-d3bmj/def2-tzvp")
#
#### matchMinima prep and run
#for adir, cdir, ddir, mmdir, smi in zip(adirs, cdirs, ddirs, mmdirs, smiles):
#    osdf = smi.split('.')[0] + '-optimized.sdf'
#    newname = ""
#    matchMinima.prep(adir, osdf, mmdir, newname)
#
#for mmdir, ref, sdflist in zip(mmdirs, refsdfs, sdflists):
#    matchMinima.matchMinima(mmdir, ref, sdflists)
