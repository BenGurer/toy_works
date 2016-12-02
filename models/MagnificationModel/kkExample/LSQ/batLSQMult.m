function [Ctip,Ctail,W,RMSD,ad2,ad3] = batLSQMult(PF,af,Lev,Al,M,y)

x = lsqnonlin(@(x)funSqEMult(PF,af,Lev,Al,M,exp(mean(log(y),2)),x),[1 0.1 -10],[0.1 0.01 -50],[10 1 0]);
Ctip = x(1); Ctail = x(2); W = x(3); 

ePProbe = funROEX(PF,Lev,x(1),x(2),x(3),false);
PR = sum(ePProbe.excInt);

ad2 = zeros(1,length(af));
ad3 = zeros(1,length(af));
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,x(1),x(2),x(3),false);
    supp2 = (1-exp(-M*1000)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al);
    supp3 = (1-exp(-M*1500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*1000)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al);
    ad2(I) = (1-sum(supp2.*ePProbe.excInt)/PR)*100;
    ad3(I) = (1-sum(supp3.*ePProbe.excInt)/PR)*100;
end
RMSD = sqrt(mean(([ad2 ad3]-[y(:,1)' y(:,2)']).^2));


