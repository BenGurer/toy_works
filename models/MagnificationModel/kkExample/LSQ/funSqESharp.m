function sqe = funSqESharp(PF,af,Lev,y,x)

ePProbe = funROEX(PF,Lev,1,x(1),x(2),false);
PR = sum(ePProbe.excInt);

sqe = zeros(1,length(af));
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,1,x(1),x(2),false);
    ePProbe = funROEXSharp(PF,Lev,1,x(1),x(2),(ePAdapt.excInt/max(ePAdapt.excInt)).^x(3),x(4),false);
    sqe(I) = (1-sum(ePProbe.excInt)/PR)*100-y(I);
end
