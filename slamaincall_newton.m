function[tvec,Pfull,HLfull,wfull,exitflagmat,itermat]=slamaincall_newton(par,tend,deltat,node_coords,...
    el_node_map,bcnodes)
%-------------------------------------------------------------------------%
% slamaincall_newton solves the system of equations 
% $\mathbf{A}\Delta d =-\mathbf{r}$
% for a given set of process parameters and geometry specifications.
% [t,Pfull,HLfull,wfull]=slamaincall(par,tend,node_coords,el_node_map,bcnodes) 
% computes Pfull (pressure p(x,t)),  HLfull (gap height h(x,t)) & wfull (film 
% deformation u_z(x,t)) for a given set of process parameters (par), time span 
% ([0,tend]) and a part geometry specified by set of nodes (node_coords) and
% elements (el_node_map) . 
% 
% INPUTS
% par : (6,1) Column vector of process parameters [V; mu; EF; NUF; HF0; h0] where 
% V denotes platform velocity, mu denotes liquid viscosity, EF denotes Youngs 
% modulus of film, NUF denotes poisson ratio of the film, HF0 denotes the initial 
% film thickness and h0 denotes liquid thickness.    
% node_coords: (N,3) matrix [i,x_i,y_i], i=1,...N 
% el_node_map: (M,5) matrix [J,j_1,j_2,j_3,j_4], J=1,..,M 
%  bcnodes : (N_b,1) column vector of  node indices of Dirichlet BC. 
%-------------------------------------------------------------------------%
nonbcnodes=node_coords(:,1);
n=length(node_coords(:,1));
nonbcnodes(bcnodes)=[];
%% FE global matrices
[Kglobal,Mglobal,bvec]=kmglobalfn(el_node_map,node_coords);
%Apply BC
Kglobal(bcnodes,:)=[];
Kglobal(:,bcnodes)=[];
Mglobal(bcnodes,:)=[];
Mglobal(:,bcnodes)=[];
bvec(bcnodes)=[];
% Process Parameters 
V=par(1);
mu=par(2);
EF=par(3);
NUF=par(4);
HF0=par(5);
h0=par(6);
[g0,g1,g2]=slaprocess(par); %$\alpha,\beta,\gamma$
gvec=[g0,g1,g2];
% Initial Conditions
ypr0=zeros(length(nonbcnodes),1);
yhl0=h0*ones(length(nonbcnodes),1);
y0=[yhl0;ypr0];
tvec=deltat:deltat:tend;
y=zeros(length(tvec),length(y0));
hlind=1:length(nonbcnodes);
pind=hlind+length(nonbcnodes);
Yold=y0;
opts=optimoptions('fsolve','Algorithm','trust-region-reflective',...
    'Jacobian','on','FunValCheck','off','MaxIter',50,'Display','off',...
    'TolX',5e-4,'TolFun',5e-4);
[GPX4,GPY4,GWEI4]=gaussint(4); % Gauss quadratue points and weights
[PHIMAT4,PHICMAT4,PHIEMAT4]=isoparshapefn(GPX4,GPY4,GWEI4); %$\mathbf{B^(4)}
[PHIXCELL,PHIYCELL,WTCELL]=knlcellfn(el_node_map,node_coords);
 for tind=1:length(tvec)
     tind
     told=tvec(tind);
 [Ynew,fval,exitflag,output] = fsolve(@(Y)slanewtfn(Y,Yold,told,deltat,par,...
     gvec,node_coords,el_node_map,Kglobal,Mglobal,bvec,bcnodes,nonbcnodes,...
     PHIMAT4,PHIXCELL,PHIYCELL,WTCELL),Yold,opts);
 y(tind,:)=Ynew;
 Yold=Ynew;
 outputcell{tind}=output;
 itermat(tind)=output.iterations;
 exitflagmat(tind)=exitflag;
 end
Pfull=zeros(length(tvec),size(node_coords,1));
HLfull=zeros(length(tvec),size(node_coords,1));
wfull=zeros(length(tvec),size(node_coords,1));
for tind=1:length(tvec)
Pfull(tind,nonbcnodes)=y(tind,pind);
HLfull(tind,bcnodes)=h0+V*tvec(tind);
HLfull(tind,nonbcnodes)=y(tind,hlind);
wfull(tind,:)=(h0+V*tvec(tind))*ones(1,n)-HLfull(tind,:);
end