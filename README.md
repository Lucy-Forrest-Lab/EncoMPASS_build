# EncoMPASS
Source code for building the EncoMPASS database

## Contributors
Antoniya A. Aleksandrova  
Edoardo Sarti  
Lucy R. Forrest  


## Release Notes (CHANGELOG) v1.1.0
Updated API calls to OPM as of 2025?  
Updated API calls to PDB as of 2025?   
Compatible with python 3.12  
Fixed bug with handling Muscle and updated to Muscle v5  
Fixed bug in output folder specification in `complete_information`   
Created wrapper code `run_encompass.py` to allow dataset compilation to be run in stages  
Updated to newer version of PPM ?  


## Release Notes v1.0.0
#### Structure Retrieval
- added information about processing and decision-making steps to the header of each structure
#### Structure Alignment
- TMs of all sequence-related chains are considered when deciding which comparisons to make
- 1 & 2 TM chains have a different set of rules from larger chains, which includes a condition about the size of the domains on either side of the membrane
#### Symmetry Algorithms Used
 - CE-Symm v2.2.3
 - QuatSymm v2.2.3
 - SymD v1.6
 - AnaNaS v1.1
#### Multi-step Symmetry Selection (MSSD)
- We've integrated QuatSymm in the MSSD procedure. We postprocess QuatSymm results to guess the specific repeat range and use the output only if the resulting symmetry has comparable RMSD and TM-score to the one reported by QuatSymm. 
- Quaternary symmetries with only 1 TM chain in a repeat are now considered acceptable and are reported.

## Dependencies & Containers
The public version container (from 2203-12-19) is available here: https://github.com/Lucy-Forrest-Lab/EncoMPASS-containers-deps
To summarize its content:
PPM v2.0 is used to insert structures in the membrane if the OPM structure for the associated biological assembly is not available
[PPM website](https://opm.phar.umich.edu/ppm_server2_cgopm)   
MUSCLE v3.8.31 is used to align sequences  
FrTMAlign is used to align structures. See [FrTM-Align website](https://sites.gatech.edu/cssb/fr-tm-align/)  
Symmetries are calculated with:
 - CE-Symm v2.2.3
 - QuatSymm v2.2.3
 - SymD v1.61 and v1.3w
 - AnaNaS v1.1

