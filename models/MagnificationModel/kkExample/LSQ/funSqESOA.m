function sqe = funSqESOA(PF,af,Lev,Ctail,W,Al,y,x)

ePProbe = funROEX(PF,Lev,1,Ctail,W,false);
PR = sum(ePProbe.excInt);

sqe = zeros(1,length(af));
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,1,Ctail,W,false);
    ePProbe = funROEXSharp(PF,Lev,1,Ctail,W,(ePAdapt.excInt/max(ePAdapt.excInt)).^Al,x,false);
    sqe(I) = (1-sum(ePProbe.excInt)/PR)*100-y(I);
end
