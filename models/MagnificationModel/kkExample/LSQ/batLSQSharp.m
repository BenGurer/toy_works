function [Ctip,Ctail,W,Al,s,RMSD,ad] = batLSQSharp(PF,af,Lev,y,y0)

x = lsqnonlin(@(x)funSqESharp(PF,af,Lev,y0,x),[0.1 -10 1 5],[0.01 -25 0 1],[1 0 1 15]);
Ctip = 1; Ctail = x(1); W = x(2); Al = x(3); 

ePProbe = funROEX(PF,Lev,1,Ctail,W,false);
PR = sum(ePProbe.excInt);

s = zeros(1,size(y,2));
ad = zeros(length(af),size(y,2));
for I = 1:size(y,2)
    s(I) = lsqnonlin(@(x)funSqESOA(PF,af,Lev,Ctail,W,Al,y(:,I),x),5,1,15);
    for II = 1:length(af)
        ePAdapt = funROEX(af(II),Lev,Ctip,Ctail,W,false);
        ePProbe = funROEXSharp(PF,Lev,Ctip,Ctail,W,(ePAdapt.excInt/max(ePAdapt.excInt)).^Al,s(I),false);
        ad(II,I) = (1-sum(ePProbe.excInt)/PR)*100;
    end
end
RMSD = sqrt(mean((ad(:)-y(:)).^2));




