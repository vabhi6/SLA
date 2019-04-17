function[PHIXMAT,PHIYMAT,WTVEC]=shapeder(GPX,GPY,GWEI,PHIMAT,PHICMAT,PHIEMAT,nodeCOORDS)
%% Inputs %%
% -------------------------------------------------------------------------
% GWEI: (9,1) vector of 2d gauss weights 
% PHIMAT: Matrix of shape functions (9,8)
%   PHIMAT(i,j)=N_j(xi(i),eta(i)),  j=1 to 8
% PHICMAT, PHIEMAT: Matrices of shape functions derivatives (9,8)
%   PHICMAT(i,j)={D(N_j)/(D xi)}(xi(i),eta(i)),i=1 to 9, j=1 to 8 
%   PHIEMAT(i,j)={D(N_j)/(D eta)}(xi(i),eta(i)),i=1 to 9, j=1 to 8 
%   PHICMAT(i,j)={D(N_j)/(D xi)}(xi(i),eta(i)),i=1 to 9, j=1 to 8 
%   PHIEMAT(i,j)={D(N_j)/(D eta)}(xi(i),eta(i)),i=1 to 9, j=1 to 8 
% nodeCoords: Matrix of nodal coordinates (8,2)
%   nodeCoords(i,:)=[xi(i),eta(i)], Matrix of material points
% -------------------------------------------------------------------------
% Initialize variables
nn=4; % Number of nodes per element
ngpt=length(GPX);
PHIXMAT=zeros(ngpt,nn);
PHIYMAT=zeros(ngpt,nn);

% Isoparametric derivatives of coordinates
DXDC=PHICMAT*nodeCOORDS(:,1);
DXDE=PHIEMAT*nodeCOORDS(:,1);
DYDC=PHICMAT*nodeCOORDS(:,2);
DYDE=PHIEMAT*nodeCOORDS(:,2);

% Element jacobian
AJACOB=(DXDC.*DYDE-DXDE.*DYDC);

for j=1:ngpt
PHIXMAT(j,:)=(PHICMAT(j,:)*DYDE(j)-PHIEMAT(j,:)*DYDC(j))/AJACOB(j);
PHIYMAT(j,:)=(PHIEMAT(j,:)*DXDC(j)-PHICMAT(j,:)*DXDE(j))/AJACOB(j);
end


WTVEC=GWEI.*AJACOB;

