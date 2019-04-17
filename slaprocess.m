function[g0,g1,g2]=slaprocess(par)
%--------------------------------------------------------------------------
% slaprocess Functions of poisson's ratio for a given set of process
% parameters (par).
%--------------------------------------------------------------------------
V=par(1);
mu=par(2);
EF=par(3);
NUF=par(4);
HF0=par(5);
h0=par(6);

f0=(1+NUF)*(1-2*NUF)/(EF*(1-NUF));
f1=NUF*(1+NUF)*(1-4*NUF)/(3*EF*(1-NUF)^2);
f2=((1+NUF)*(1-4*NUF)/(2*EF*(1-NUF)));

g0=f0*HF0;
g1=-(f1)*HF0^3;
g2=(-f2*HF0^2)/(2);

