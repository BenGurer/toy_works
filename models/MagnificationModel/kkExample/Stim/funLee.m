function lee = funLee(F,SF,Nfft)

f = SF/2*linspace(0,1,Nfft); DF = diff(f(1:2));

ee = 1./lcfErb(f);

F1 = lcfInvNErb(lcfNErb(F)-0.5);
F2 = lcfInvNErb(lcfNErb(F)+0.5);
bp = zeros(1,Nfft);
bp(round(F1/DF):round(F2/DF)) = 1;

lee(1) = 10*log10(sum(ee.*bp)/sum(ee));
lee(2) = 10*log10(lcfIntErb(F1,F2)/lcfIntErb(0,SF/2));

% ***** lcfErb *****
function erb = lcfErb(f)

A = 24.7/1000; B = 4.37;
erb = A*(B*f+1);

% ***** lcfNErb *****
function nerb = lcfNErb(f)

A = 24.7/1000; B = 4.37;
nerb = 1/(A*B)*log(B*f+1);

% ***** lcfInvNErb *****
function f = lcfInvNErb(nerb)

A = 24.7/1000; B = 4.37;
f = 1/B*(exp(A*B*nerb)-1);

% ***** lcfIntErb *****
function I = lcfIntErb(Fl,Fu)

A = 24.7/1000; B = 4.37;
I = 1/(A*B)*(log(B*Fu+1)-log(B*Fl+1));



