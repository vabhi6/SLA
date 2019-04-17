function [Kglobal,Mglobal,bvec]=kmglobalfn(el_node_map,node_coords)
%--------------------------------------------------------------------------
%kmglobalfn computes the matrices(K,M) and column vector(f) for a
%prescribed finite element mesh. 
%--------------------------------------------------------------------------

node_vec=node_coords(:,1); % Vector of node labels
n=length(node_vec); % Total Number of nodes
nel=length(el_node_map(:,1));
nn=4;% number of nodes per element
ntriplets = nel*nn^2 ;
I = zeros (ntriplets, 1) ;
J = zeros (ntriplets, 1) ;
KX = zeros (ntriplets, 1) ;
MX = zeros (ntriplets, 1) ;
bvec=zeros (n, 1) ;
[GPX,GPY,GWEI]=gaussint(2);
ngpt=length(GPX);
[PHIMAT,PHICMAT,PHIEMAT]=isoparshapefn(GPX,GPY,GWEI);
 ntriplets=0;
for i=1:nel
    nind=el_node_map(i,2:5);
    elnodeCOORDS=node_coords(nind,2:3);
    [PHIXMAT,PHIYMAT,WTVEC]=shapeder(GPX,GPY,GWEI,PHIMAT,PHICMAT,...
        PHIEMAT,elnodeCOORDS);
    Kel=zeros(nn);
    Mel=zeros(nn);
    Bel=zeros(nn,1);
    for j=1:ngpt
        Kel=Kel+WTVEC(j)*(PHIXMAT(j,:)'*PHIXMAT(j,:) + PHIYMAT(j,:)'*PHIYMAT(j,:));
        Mel=Mel+WTVEC(j)*PHIMAT(j,:)'*PHIMAT(j,:);
        Bel=Bel+WTVEC(j)*PHIMAT(j,:)';
    end

    belvec=zeros(n,1);
    belvec(nind)=Bel;
    bvec=bvec+belvec;
    for krow=1:nn
            for kcol=1:nn
                ntriplets = ntriplets + 1 ;
                I (ntriplets) = nind (krow) ;
                J (ntriplets) = nind (kcol) ;
                KX (ntriplets) = Kel (krow,kcol) ;
                MX (ntriplets) = Mel (krow,kcol) ;
            end
    end
end
Kglobal = sparse (I,J,KX,n,n);
Mglobal = sparse (I,J,MX,n,n);