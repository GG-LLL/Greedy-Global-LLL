# Greedy-Global-LLL

This repository contains implementations accompanying the paper [5] titled "A Greedy Global Framework for LLL".
 The algorithms implemented are: Pot-DeepLLL [2], SS-DeepLLL [3], and their Greedy-Global variants proposed in [5].
 The Pot-DeepLLL, SS-DeepLLL, Pot-GGLLL and SS-GGLLL algorithms use the GSO update techniques from [4], that make them more efficient than BKZ-8!
 SS-GGLLL gives shorter vectors than BKZ-12 for dimensions 100 and after.

The descriptions of the files that can be found in this repository are listed below.

- common.cpp &ensp; common.h &ensp; standalone_reduction.cpp &ensp; preprocessed_reduction.cpp
For a given input basis, the functions contained in common.cpp may be used to run various lattice basis reduction algorithms. The menu-driven execution of individual algorithms is in standalone_reduction.cpp and preprocessed_reduction.cpp for multi-precision and standard datatype implementations respectively.
  - LLL (multi-precision implementation with overestimated precision)
  - Pot-DeepLLL (standard datatype implementation and multi-precision implementation with overestimated precision)
  - SS-DeepLLL (standard datatype implementation and multi-precision implementation with overestimated precision)
  - Pot-GGLLL (standard datatype implementation and multi-precision implementation with overestimated precision)
  - SS-GGLLL (standard datatype implementation and multi-precision implementation with overestimated precision)

- generate_SVP_bases.cpp  
Can be used to generate SVP Challenge style bases. 

- collate.cpp  
Can be used to collate the results from multiple output files in one place.

- SVP_Bases &ensp; SVP_Bases_FPLLL &ensp; SVP_Bases_LLLReduced   
The input bases are contained in the zip file SVP_Bases. The same bases are in the zip file SVP_Bases_FPLLL, formatted for the use of FPLLL functions. The zip file SVP_Bases_LLLReduced contains the SVP_Bases in the other files, but 0.99-LLL reduced.

- SVP_Outputs_Standalone &ensp; SVP_Outputs_Preprocessed  
The outputs from our programs are contained in the zip files SVP_Outputs_Standalone and SVP_Outputs_Preprocessed for multi-precision and standard implementations respectively.

- Standalone_Averages &ensp; Preprocessed_Averages  
This folder contains the summaries reported in our paper.

References:

 [1] Claus-Peter Schnorr and Martin Euchner. Lattice basis reduction: improved practical algorithms and solving subset sum problems. Mathematical programming, 66(1):181–199, 1994.
 
 [2] Felix Fontein, Michael Schneider, and Urs Wagner. PotLLL: a polynomial time version of LLL with deep insertions. Designs, Codes and Cryptography, 73(2):355–368, 2014.
 
 [3] Masaya Yasuda and Junpei Yamaguchi. A new polynomial-time variant of LLL with deep insertions for decreasing the squared-sum of Gram–Schmidt lengths. Designs, Codes and Cryptography, 87(11):2489–2505, 2019.
 
 [4] Junpei Yamaguchi and Masaya Yasuda. Explicit formula for gram-schmidt vectors in LLL with deep insertions and its applications. In Jerzy Kaczorowski, Josef Pieprzyk, and Jacek Pomykala, editors, Number-Theoretic Methods in Cryptology, pages 142–160, Cham, 2018. Springer International Publishing.
 
 [5] Sanjay Bhattacherjee, Julio Hernandez-Castro and Jack Moyler. A Greedy Global Framework for LLL, Cryptology ePrint Archive, Paper 2023/261, https://eprint.iacr.org/2023/261.

Contributors:
 Sanjay Bhattacherjee (sanjay.bhattacherjee@gmail.com) and Jack Moyler (j.d.moyler@gmail.com)
