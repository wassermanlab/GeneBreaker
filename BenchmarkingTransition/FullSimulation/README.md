# Full Simulation

> The purpose of this pipeline is to simulate a full genome fastq file, with the desired variants incorporated.

## Overview

1. Acquire VarSim, and some extra data files:

```
sh DeployVarsimAndTest.sh 
```

2. Run a simulation with example variants, but make sure you change the paths within this script:

```
sh VarSim_GeneBreakerLargeChallengingNoncodingMEIs.sh
```

The output of this will include: 
- A directory from VarSim with the raw fastq data, and details of the simulation.
- A mapped BAM file of your aligned reads. (NOTE: If you don't want this, just remove from the script).


3. (Optional) Create your own example run by inserting your own VCFs into the VarSim parameter '--vcfs'




