#!/usr/bin/env python

import os, sys
import smi2confs
import filterConfs
import confs2psi
import getPsiResults
import matchMinima

adirs = []
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-A' ]  
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-A'  ] 
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-A']  
adirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-A']  

cdirs = []
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-C' ]  
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-C'  ] 
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-C']  
cdirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-C']  

ddirs = []
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/1_chains-D' ]  
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/2_rings-D'  ] 
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/3_sterics-D']  
ddirs += ['/work/cluster/limvt/qm_AlkEthOH/pipeline/4_diverse-D']  

smiles = []
smiles += ['AlkEthOH_chain_tiny.smi' ] 
smiles += ['AlkEthOH_rings_tiny.smi' ] 
smiles += ['AlkEthOH_rings_supp1.smi'] 
smiles += ['diverse.smi'             ] 


#### For the A directories, MM opt & filter
#for wdir, smi in zip(adirs, smiles):
#    print "\n", wdir, smi
#    sdf = smi.split('.')[0] + '.sdf'
#    smi2confs.smi2confs(wdir, smi)
#    filterConfs.filterConfs(wdir, sdf, "MM Szybki SD Energy")

#### For the C directories, no MM opt & filter
#for wdir, smi in zip(cdirs, smiles):
#    print "\n", wdir, smi
#    sdf = smi.split('.')[0] + '.sdf'
#    smi2confs.smi2confs(wdir, smi, True, False)
#    filterConfs.filterConfs(wdir, sdf, "MM Szybki Single Point Energy")

### For A and C directories, confs2psi. Then getPsiResults, filter (after QM).
for i, smi in enumerate(smiles):
    print "\n", i, smi
    sdf = smi.split('.')[0] + '-minima.sdf'
    adir = adirs[i]
    cdir = cdirs[i]
    confs2psi.confs2psi(adir, sdf, 'mp2','def2-sv(p)', False, "1.5 Gb") 
#    confs2psi.confs2psi(cdir, sdf, 'mp2','def2-sv(p)', False, "1.5 Gb") 


#### For the D directories, confs2psi with green star method.
#for wdir, smi in zip(ddirs, smiles):
#    sdf = smi.split('.')[0] + '-minima.sdf'
#    confs2psi.confs2psi(wdir, sdf, 'b3lyp-d3mbj','def2-tzvp', False, "1.5 Gb") 
