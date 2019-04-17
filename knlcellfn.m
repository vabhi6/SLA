function [PHIXCELL,PHIYCELL,WTCELL]=knlcellfn(el_node_map,node_coords)
%--------------------------------------------------------------------------
%knlcellfn computes the derivatives of shape functions cooresponding to each element 
%$\mathbf{B}^(4)$ for a prescribed finite element mesh. 
%--------------------------------------------------------------------------
node_vec=node_coords(:,1); % Vector of node labels
n=length(node_vec); % Total Number of nodes
nel=length(el_node_map(:,1));
nn=4;% number of nodes per element
[GPX,GPY,GWEI]=gaussint(4);
ngpt=length(GPX);
[PHIMAT,PHICMAT,PHIEMAT]=isoparshapefn(GPX,GPY,GWEI);

for i=1:nel
    nind=el_node_map(i,2:5);
    elnodeCOORDS=node_coords(nind,2:3);
    [PHIXMAT,PHIYMAT,WTVEC]=shapeder(GPX,GPY,GWEI,PHIMAT,PHICMAT,...
        PHIEMAT,elnodeCOORDS);  
PHIXCELL{i}=PHIXMAT;
PHIYCELL{i}=PHIYMAT;
WTCELL{i}=WTVEC;
end

    

