# SLA
Finite Element code to compute separation force on parts printed using constrained-surface stereolithography (SLA)
Problem: To compute lubrication forces on a rigid cylinder separating from a thin soft planar coating (bonded to a rigid substrate) in a liquid medium
Main script: slafemmain.m
Inputs Files: Specify three text files (available in a suitable path) containing:
(i) nodal indices and their coordinates  
(ii) Element connectivity matrix (quadilateral elements)
(iii) Dirichlet boundary condition nodal indices
Input parameters:
v_z % Separation rate
mu_l % Liquid viscosity
E_f % Youngs Modulus of film coating 
nu_f % Poisson's ratio of film coating
h_f0 % Film thickness
h_l0 % Initial liquid height (assumed to be uniform spatially)
Output variables:
(i) Nodal pressures at time instants
(ii) Nodal heights at time instants
