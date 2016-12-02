function [Ctip,Ctail,W,Al,RMSD,ad] = batLSQSingle(PF,af,Lev,soa,M,y,y0)

x = lsqnonlin(@(x)funSqESingle(PF,af,Lev,y0,x),[0.1 -10 1],[0.01 -25 0],[1 0 1]);
Ctip = 1; Ctail = x(1); W = x(2); Al = x(3);

ePProbe = funROEX(PF,Lev,Ctip,Ctail,W,false); PR = sum(ePProbe.excInt); 
ad = zeros(length(af),length(soa));
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,Ctip,Ctail,W,false); 
    for II = 1:length(soa)
        supp = 1-exp(-M*soa(II))*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al; 
        ad(I,II) = (1-sum(supp.*ePProbe.excInt)/PR)*100; 
    end
end
RMSD = sqrt(mean((ad(:)-y(:)).^2));

    
