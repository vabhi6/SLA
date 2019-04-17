function[F,Jacmat]=slanewtfn(Y,Yold,told,deltat,par,gvec,node_coords,...
    el_node_map,Kglobal,Mglobal,bvec,bcnodes,nonbcnodes,PHIMAT4,...
    PHIXCELL,PHIYCELL,WTCELL)
%--------------------------------------------------------------------------
% slanewtfn computes the residual $\mathbf{r}$ and Jacobian $\mathbf{A}$
%--------------------------------------------------------------------------
%t denotes t+dt
t=told+deltat;
V=par(1);
h0=par(6);
phi=h0+V*t;
beta=phi*gvec(3);
n_nbc=length(nonbcnodes);
Ydot=(Y-Yold)/deltat;
hdot=Ydot(1:n_nbc);
pdot=Ydot(n_nbc+1:2*n_nbc);
[Knlglobal,Anlglobal]=Knlcalc2(t,Y,par,gvec,node_coords,el_node_map,...
    Kglobal,Mglobal,bvec,bcnodes,nonbcnodes,PHIMAT4,...
    PHIXCELL,PHIYCELL,WTCELL);
Cglobal=-gvec(1)*Mglobal-gvec(2)*Kglobal;
F=[zeros(n_nbc) Knlglobal;Mglobal Cglobal]*Y+[Mglobal*hdot+ Anlglobal*pdot;-phi*bvec];
Jacmat=[Mglobal/deltat Anlglobal/deltat+Knlglobal;Mglobal Cglobal];
