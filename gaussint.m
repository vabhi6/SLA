function[GPX,GPY,GWEI]=gaussint(order)
if(order==1)
    gauss_pt=0;
    GWE=2;
    GPX=gauss_pt;
    GPY=gauss_pt;
    GWEI=GWE*GWE;
elseif(order==2)
    gauss_pt=[-1,1]*sqrt(1/3);
    GPX=zeros(4,1);
    GPY=zeros(4,1);
    count=0;
    for i=1:2
        for j=1:2
            count=count+1;
            GPX(count)=gauss_pt(i);
            GPY(count)=gauss_pt(j);
        end
    end
    GWE=[1,1];
    GWEI= reshape(GWE'*GWE,4,1);
    
elseif(order==3)
    gauss_pt=[-1,0,1]*sqrt(6/10);
    GPX=zeros(9,1);
    GPY=zeros(9,1);
    count=0;
    for i=1:3
        for j=1:3
            count=count+1;
            GPX(count)=gauss_pt(i);
            GPY(count)=gauss_pt(j);
        end
    end
    GWE=[5/9, 8/9, 5/9];
    GWEI= reshape(GWE'*GWE,9,1);
    
elseif(order==4)
    r=sqrt(1.2);
    gauss_pt1=sqrt((3+2*r)/7);
    gauss_pt2=sqrt((3-2*r)/7);
    gauss_pt=[-gauss_pt1,gauss_pt1,-gauss_pt2,gauss_pt2];
    GPX=zeros(16,1);
    GPY=zeros(16,1);
    count=0;
    for i=1:4
        for j=1:4
            count=count+1;
            GPX(count)=gauss_pt(i);
            GPY(count)=gauss_pt(j);
        end
    end
    GWE=[0.5-(1/(6*r)),0.5-(1/(6*r)),0.5+(1/(6*r)),0.5+(1/(6*r))];
    GWEI= reshape(GWE'*GWE,16,1);
end
