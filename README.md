# EncoMPASS
Source code for building the EncoMPASS database

## Contributors
Antoniya A. Aleksandrova  
Edoardo Sarti  
Lucy R. Forrest  

## Release Notes & Current Updates
### Current Updates
#### Structure Retrieval
- added information about processing and decision-making steps to the header of each structure
#### Structure Alignment
- TMs of all sequence-related chains are considered when deciding which comparisons to make
- 1 & 2 TM chains have a different set of rules from larger chains, which includes a condition about the size of the domains on either side of the membrane
#### Symmetry Algorithms Used
 - CE-Symm v2.2.2
 - QuatSymm v2.2.2
 - SymD v1.6
 - AnaNaS v1.1  (not currently used for analysis, but available)

#### Multi-step Symmetry Selection (MSSD)
- We've integrated QuatSymm in the MSSD procedure. We postprocess QuatSymm results to guess the specific repeat range and use the output only if the resulting symmetry has comparable RMSD and TM-score to the one reported by QuatSymm. 
- Quaternary symmetries with only 1 TM chain in a repeat are now considered acceptable and are reported.


## Dependencies & Containers
PPM v2.0 is used to insert structures in the membrane if the OPM structure for the associated biological assembly is not available. For more information, please see the
[PPM website](https://opm.phar.umich.edu/ppm_server2_cgopm)   
MUSCLE v3.8.3 is used to align sequences. See [Muscle v3 website](https://drive5.com/muscle/downloads_v3.htm)  
FrTMAlign is used to align structures. See [FrTM-Align website](https://sites.gatech.edu/cssb/fr-tm-align/)  
Symmetries are calculated with: 
 - [CE-Symm v2.2.2](https://github.com/rcsb/symmetry/releases/tag/symmetry-2.2.2)
 - [QuatSymm v2.2.2](https://github.com/rcsb/symmetry/releases/tag/symmetry-2.2.2)
 - SymD v1.6
 - [AnaNaS v1.1](https://team.inria.fr/nano-d/software/ananas)

These dependencies can be copied to and run from a container, instructions for which are available on [github](https://github.com/Lucy-Forrest-Lab/EncoMPASS-containers-deps)


 


