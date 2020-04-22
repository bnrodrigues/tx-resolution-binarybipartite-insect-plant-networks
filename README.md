# tx-resolution-binarybipartite-insect-plant-networks
Rodrigues, B.N. & Boscolo, D. (2020) Do bipartite binary antagonistic and mutualistic networks have different responses to the taxonomic resolution of nodes? Ecological Entomology. DOI:10.1111/een.12844

In order to run the first step ("_step1_several-webs_metrics.R"), you need to have the set of networks (sources decribed in Table S1 and Table S2), with insects in the columns and plants in the rows. You need to put the antagosnic matrixes in the folder "antagonism-all" and the mutulistic ones in "mutualism-all". Both folders need to be in the folder "adjMatrices". You also can change the script to work with your the directory strucuture.
Our folder structure was: 
├───adjMatrices
│   ├───antagonism-all
│   └───mutualism-all
├───output
└───scriptR (this was our working directory)

The files without the prefix "_step" are graphic and table of the results generated from the steps 1 to 3.
