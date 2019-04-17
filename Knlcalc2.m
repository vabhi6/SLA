function[Knlglobal,Anlglobal]=Knlcalc2(t,Y0,par,gvec,node_coords,...
    el_node_map,Kglobal,Mglobal,bvec,bcnodes,nonbcnodes,PHIMAT4,...
    PHIXCELL,PHIYCELL,WTCELL)
%-------------------------------------------------------------
% Knlcalc2 computes the matrices $\mathbf{G}$ and $\mathbf{D}$
%-------------------------------------------------------------

tnew=t;
V=par(1);
mu=par(2);
EF=par(3);
NUF=par(4);
HF0=par(5);
h0=par(6);

%% Assemble K_{nl}

n=size(node_coords,1); % Total Number of nodes

nel=size(el_node_map,1);
nn=4;% Number of nodes per element

ntriplets = nel*nn^2;

phi=h0+V*tnew;

Y0_hl=Y0(1:length(nonbcnodes));
Hlfull=phi*ones(n,1);
Hlfull(nonbcnodes)=Y0_hl;
betasc=phi*gvec(3);
I = zeros (ntriplets, 1) ;
J = zeros (ntriplets, 1) ;
KnlX = zeros (ntriplets, 1);
AnlX = zeros (ntriplets, 1);

ntriplets=0;
for i=1:nel
    
    nind=el_node_map(i,2:5);
    elnodeCOORDS=node_coords(nind,2:3);
    PHIXMAT=PHIXCELL{i};
      PHIYMAT=PHIYCELL{i};
        WTVEC=WTCELL{i};
    Hlel=Hlfull(nind);
    HLvec=PHIMAT4*Hlel;
    alpha=HLvec.^3/(12*mu);
    beta=betasc*HLvec;
      if(sum(alpha<0)>0)
      disp(i);
      end
    Knlel=zeros(nn);
    Anlel=zeros(nn);
        for j=1:length(WTVEC)
        Knlel=Knlel+WTVEC(j)*alpha(j)*(PHIXMAT(j,:)'*PHIXMAT(j,:) + PHIYMAT(j,:)'*PHIYMAT(j,:));
        Anlel=Anlel+WTVEC(j)*beta(j)*(PHIXMAT(j,:)'*PHIXMAT(j,:) + PHIYMAT(j,:)'*PHIYMAT(j,:));
        end
         for krow=1:nn
            for kcol=1:nn
                ntriplets = ntriplets + 1 ;
                I (ntriplets) = nind (krow) ;
                J (ntriplets) = nind (kcol) ;
                KnlX (ntriplets) = Knlel (krow,kcol) ;
                AnlX (ntriplets) = Anlel (krow,kcol) ;
            end
         end
end
Knlglobal = sparse (I,J,KnlX,n,n) ;
Anlglobal = sparse (I,J,AnlX,n,n) ;
Knlglobal(bcnodes,:)=[];
Knlglobal(:,bcnodes)=[];
Anlglobal(bcnodes,:)=[];
Anlglobal(:,bcnodes)=[];
