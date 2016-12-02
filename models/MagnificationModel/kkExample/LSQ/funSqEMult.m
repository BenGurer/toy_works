function sqe = funSqEMult(PF,af,Lev,Al,M,y,x)

ePProbe = funROEX(PF,Lev,x(1),x(2),x(3),false);
PR = sum(ePProbe.excInt);

sqe = zeros(1,length(af));
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,x(1),x(2),x(3),false);
    supp2 = (1-exp(-M*1000)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al);
    supp3 = (1-exp(-M*1500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*1000)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al);
    Ad2 = (1-sum(supp2.*ePProbe.excInt)/PR)*100;
    Ad3 = (1-sum(supp3.*ePProbe.excInt)/PR)*100;
    sqe(I) = exp(mean(log([Ad2 Ad3])))-y(I);
end
