function[PHIMAT,PHICMAT,PHIEMAT]=isoparshapefn(GPX,GPY,GWEI)
ngpt=length(GPX);
nn=4;
PHIMAT=zeros(ngpt,nn);
PHICMAT=zeros(ngpt,nn);
PHIEMAT=zeros(ngpt,nn);

for i=1:ngpt 
    C=GPX(i);
    E=GPY(i);
    N1=0.25*(1-C)*(1-E);    
    N2=0.25*(1+C)*(1-E);
    N3=0.25*(1+C)*(1+E);
    N4=0.25*(1-C)*(1+E);
    PHIMAT(i,:)=[N1,N2,N3,N4];
    PHICMAT(i,:)=[ E/4 - 1/4, 1/4 - E/4, E/4 + 1/4, - E/4 - 1/4];
    PHIEMAT(i,:) = [ C/4 - 1/4, - C/4 - 1/4, C/4 + 1/4, 1/4 - C/4];
end