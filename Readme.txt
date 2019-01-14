TVL0_select

*****************************************************************************************************************
* author: Jordan Frécon  											*
* institution: Univ Lyon, Ens de Lyon, Univ Claude Bernard, CNRS, Laboratoire de Physique, F-69342 Lyon, France *
* date: October 20 2016     	              									*
* License CeCILL-B                                    								*
*****************************************************************************************************************


*********************************************************
* RECOMMENDATIONS:                                   	*
* This toolbox is designed to work with Matlab 2015.a   *
*********************************************************

------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides a parameter-free denoising procedure of piecewise constant signals corrupted by Gaussian or Laplacian noise.
The corresponding method relies on an operational strategy that combines hierarchical Bayesian and TVℓ0 (also kown as ℓ2-Potts model) formulations, with the double aim of automatically tuning the regularization parameter and of maintaining computational efficiency.

This toolbox makes use of two toolboxes :
1) « Pottslab » (required) to solve the TVℓ0 minimization problem for a given regularization parameter. 
2) « Randraw » (optional) to generate Laplacian noise.
For convenience purposes, copies of the latter toolboxes are provided in the subfolders « Pottslab0.42 » and « Randraw ».

------------------------------------------------------------------------------------------------------------------------
SPECIFICATIONS for using TVL0_select:

Two demo files ‘demo_TVL0_select_Gaussian.m’ and ‘demo_TVL0_select_Laplacian.m’ are proposed depending on the Gaussian or Laplacian nature of the noise.

/!\ Run the command ‘installPottslab’ once per MATLAB session /!\

The main function is ‘TVL0_select’.

------------------------------------------------------------------------------------------------------------------------
RELATED PUBLICATION:

# J. Frecon, N. Pustelnik, N. Dobigeon, H. Wendt, and P. Abry
Bayesian selection for the l2-Potts model regularization parameter: 1D piecewise constant signal denoisingPreprint arXiv:1608.07739

------------------------------------------------------------------------------------------------------------------------