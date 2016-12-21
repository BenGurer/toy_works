function [dCtrd,dSprd,eCtrd,eSprd] = model1

N1 = lcfF2NErb(0.25); N2 = lcfF2NErb(8);
CF = 2; NOI = 5;

NRep = 100; NStim = 2.^(2:6); C = 10.^(0.25:0.25:1.25); 
NX = length(C); NY = length(NStim);

figure(1), clf, hold on
eCtrd = zeros(NY,NX); eSprd = zeros(NY,NX);
dCtrd = zeros(NY,NX); dSprd = zeros(NY,NX);
for I = 1:NY
    stimF = lcfNErb2F(N1:((N2-N1)/(NStim(I)-1)):N2);
    for II = 1:NX
        [Ctrd,Sprd,ERB1,ERB2,f,w] = lcfW(CF,C(II));
        jwd1 = zeros(1,NRep); jwd2 = zeros(1,NRep);     
        jwd3 = zeros(1,NRep); jwd4 = zeros(1,NRep);
        estERB = zeros(1,NRep);
        r = zeros(1,NStim(I));
        for III = 1:NRep
            [jwd1(III),jwd2(III),oneR] = lcfVoxResp(CF,C(II),NOI,stimF);
            jwd3(III) = abs(jwd1(III)-Ctrd)/Ctrd*100;    
            jwd4(III) = abs(jwd2(III)-Sprd)/Sprd*100;
            r = r+oneR/NRep;
        end
        eERB(I,II) = nanmean(estERB);
        eCtrd(I,II) = nanmean(jwd1); eSprd(I,II) = nanmean(jwd2);
        dCtrd(I,II) = nanmean(jwd3); dSprd(I,II) = nanmean(jwd4);
        subplot(NY,NX,(NY-I)*NX+II), hold on
        plot(f,w,'k-'), plot(stimF,r,'ro--')
        text(max(xlim),max(ylim),sprintf('Ctrd = %.2f (%.2f) kHz\nSprd = %.2f (%.2f) kHz\nERB1 = %.2f, ERB2 = %.2f',Ctrd,eCtrd(I,II),Sprd,eSprd(I,II),ERB1,ERB2),'HorizontalALignment','right','VerticalAlignment','top')
    end
end

figure(2), clf
colormap('gray')
subplot(1,2,1), contourf(C,NStim,dCtrd), set(gca,'XScale','log','YScale','log'), title('Ctrd'), xlabel('C'), ylabel('NStim')
subplot(1,2,2), contourf(C,NStim,dSprd), set(gca,'XScale','log','YScale','log'), title('Sprd'), xlabel('C')

% Local functions
function nErb = lcfF2NErb(f)
nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);

function f = lcfNErb2F(nErb)
f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;

function erb = lcfErb(f)
erb = 24.67*(4.37*f+1)/1000;

function [ECtrd,ESprd,r] = lcfVoxResp(CF,C,NOI,f)
ERB = C*lcfErb(CF);
P = 4*CF/ERB;  
g = abs(f-CF)/CF;
r = (1+P*g).*exp(-P*g)+NOI/100*randn(size(f));
ECtrd = sum(max(r,0).*f)/sum(max(r,0));
ESprd = sqrt(sum(max(r,0).*(f-ECtrd).^2)/sum(max(r,0)));
% eERB = sum(r/max(r))*diff(f(1:end));

function [Ctrd,Sprd,ERB1,ERB2,f,w] = lcfW(CF,C)
ERB1 = C*lcfErb(CF);
P = 4*CF/ERB1; 
f = CF+(-3*CF:(6*CF)/(100-1):3*CF); f = f(f>=0);
g = abs(f-CF)/CF;
w = (1+P*g).*exp(-P*g);

Ctrd = sum(w.*f)/sum(w);
Sprd = sqrt(sum(w.*(f-Ctrd).^2)/sum(w));
ERB2 = sum(w/max(w))*diff(f(1:2));




