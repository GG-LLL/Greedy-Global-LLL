# Greedy-Global-LLL

This repository contains implementations accompanying the paper [5] titled "A Greedy Global Framework for LLL".
 The algorithms implemented are: Pot-DeepLLL[2], SS-DeepLLL[3], and their Greedy-Global variants proposed in the above paper.
 The Pot-DeepLLL and SS-DeepLLL algorithms use the GSO update techniques from [4].

The descriptions of the files that can be found in this repository are listed below.

- common.cpp &ensp; common.h &ensp; individual_basis_reduction.cpp  
For a given input basis, the functions contained here may be used to run one of the 9 basis reduction algorithms listed in our paper.

- generate_SVP_basis.cpp  
Can be used to generate SVP Challenge style bases. 

- collate.cpp  
Can be used to collate the results from multiple output files in one place.

- SVP_Bases &ensp; SVP_Outputs  
The input bases are contained in the zip file SVP_Bases, while their respective outputs are contained within SVP_Outputs.

- Averages  
This folder contains the summaries reported in our paper.

References:

 [1] Claus-Peter Schnorr and Martin Euchner. Lattice basis reduction: improved practical algorithms and solving subset sum problems. Mathematical programming, 66(1):181–199, 1994.
 
 [2] Felix Fontein, Michael Schneider, and Urs Wagner. PotLLL: a polynomial time version of LLL with deep insertions. Designs, Codes and Cryptography, 73(2):355–368, 2014.
 
 [3] Masaya Yasuda and Junpei Yamaguchi. A new polynomial-time variant of LLL with deep insertions for decreasing the squared-sum of Gram–Schmidt lengths. Designs, Codes and Cryptography, 87(11):2489–2505, 2019.
 
 [4] Junpei Yamaguchi and Masaya Yasuda. Explicit formula for gram-schmidt vectors in LLL with deep insertions and its applications. In Jerzy Kaczorowski, Josef Pieprzyk, and Jacek Pomykala, editors, Number-Theoretic Methods in Cryptology, pages 142–160, Cham, 2018. Springer International Publishing.
 
 [5] Sanjay Bhattacherjee, Julio Hernandez-Castro and Jack Moyler. A Greedy Global Framework for LLL, Cryptology ePrint Archive, Paper 2023/261.

Contributors:
 Sanjay Bhattacherjee and Jack Moyler
