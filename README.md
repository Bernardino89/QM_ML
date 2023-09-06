# QM - ML
Supervised Learning for obtaining a Basis Set

****Development and supervised learning of a chemical model (Python/ Gaussian) to correct basis set effects on quantum mechanical predictions, speed up computations, and improve thermochemical accuracy.**** 

For reference, please consult the following publications:

Tirri, Bernardino, et al. "Computation of covalent and noncovalent structural parameters at low computational cost: Efficiency of the DH‐SVPD method." International Journal of Quantum Chemistry 120.13 (2020): e26233.

Li, Hanwei, et al. "Beyond chemical accuracy for alkane thermochemistry: the DH thermo approach." The Journal of Organic Chemistry 86.8 (2021): 5538-5545.

The NH3 folder pertains to the system under examination. Enclosed within this directory are all the associated files:

***optimizer.py:***

  This program has been designed to interface with the Gaussian software, enabling the development of DH-SVPD (Double-ζ) basis sets specifically. DH-SVPD leverages an error compensation     strategy to address the trade-off between basis set incompleteness errors and overlap errors (BSIE and BSSE, respectively). DH-SVPD is derived based on the zeroth-order energy exchange optimization criterion proposed by Varandas, referred to here as the constraint:
	
           _                           _  2
          |   (ED - EM)  -  (EJ + EX)   |
     min  |  _________________________  |
    {αi}  |_  (ED - EM)  +  (EJ + EX)  _|


Where the quantity E - E0 refers to the interaction energy between monomers in a molecular aggregate, and J and K represent the corresponding Coulomb and exchange energies of the dimer itself. The optimization of the exponents of the more diffuse Gaussian functions in the def2-SVPD atomic basis set, denoted here as {αi}, enables the expression of the molecular aggregate's interaction energy within its perturbative one-electron exchange development. Additionally, the program facilitates the **automated extraction of features**: the electronic energy of the dimer and monomers, MP2 energies of dimers and monomers, and exchange and correlation energies.

***starter.com:***

  Input file that encompasses the geometry, theory level, basis set for optimization, and necessary Gaussian keywords.

****slurm-123882.out:***

  Output file that contains all information.
