function[fvec]=forcecalc(t,yfull,node_coords,el_node_map)
[GPX,GPY,GWEI]=gaussint(2);
nel=length(el_node_map(:,1));
node_vec=node_coords(:,1); % Vector of node labels
n=length(node_vec); % Total Number of nodes
[PHIMAT,PHICMAT,PHIEMAT]=isoparshapefn(GPX,GPY,GWEI);
fvec=zeros(length(t),1);
for tind=1:length(t)
    force=0;
    for i=1:nel
    nind=el_node_map(i,2:5);
    elnodeCOORDS=node_coords(nind,2:3);
    elforce=elforcecalc(elnodeCOORDS,yfull(tind,nind)',GWEI,PHIMAT,PHICMAT,PHIEMAT);
    force=force+elforce;
    end
    fvec(tind)=force;
end
