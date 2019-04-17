clear all; close all;clc
% Model Inputs and Simulation parameters
v_z=0.1; 
mu_l=10*1e-8;
E_f=1;
nu_f=0.45;
h_f0=1;
h_l0=0.025;
par=[v_z;mu_l;E_f;nu_f;h_f0;h_l0]; %Model Input vector
tend=1; % Simulation Parameter $T$
deltat=0.01; % Simulation Parameter $\Delta t$
tstart=tic;
% Text files containing mesh information
node_file='node_coords.dat'; % Nodal coordinates file
el_map_file='el_node_conn.dat'; % Element Connectivity file  
bc_file='bc_nodes.dat'; % Boundary nodal indices file
% Mesh information ---------------------------------------
node_coords=importdata(node_file); % Nodal coordinates
el_node_map=importdata(el_map_file); % Element Connectivity Matrix
bcnodemat=importdata(bc_file);
bcnodevec=reshape(bcnodemat,size(bcnodemat,1)*size(bcnodemat,2),1);
bcnodevec(isnan(bcnodevec))=[];
bcnodes=sort(bcnodevec); % Dirichlet BC nodal indices
nonbcnodes=node_coords(:,1);
nonbcnodes(bcnodes)=[];
% Solve for nodal heights, pressures and out of plane film deformations 
[t,Pfull,hlfull,wfull,exitflagmat,itermat]=slamaincall_newton(par,tend,deltat,node_coords,el_node_map,bcnodes);
Pfullcell=Pfull;
wfullcell=wfull;
HLfullcell=hlfull;
[fvec]=forcecalc(t,Pfull,node_coords,el_node_map); % Force computed over a quadrant
Forcevec=-4*fvec;
telaspsed=toc(tstart);
plot(t',Forcevec)