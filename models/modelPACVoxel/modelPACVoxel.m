function [dCtrd,dSprd,eCtrd,eSprd] = modelPACVoxel (NOI)

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
        [Ctrd,Sprd,f,w] = lcfW(CF,C(II));
        jwd1 = zeros(1,NRep); jwd2 = zeros(1,NRep);     
        jwd3 = zeros(1,NRep); jwd4 = zeros(1,NRep);     
        r = zeros(1,NStim(I));
        for III = 1:NRep
            [jwd1(III),jwd2(III),oneR] = lcfVoxResp(CF,C(II),NOI,stimF);
            jwd3(III) = abs(jwd1(III)-Ctrd)/Ctrd*100;    
            jwd4(III) = abs(jwd2(III)-Sprd)/Sprd*100;
            r = r+oneR/NRep;
            fun = @lcfROEX;
xdata = stimF;
x0 = [1, 0.5];
ydata = r;
lb = [0.02, 0.1];
ub= [20, 5];
options = optimoptions('lsqcurvefit','TolFun',1e-18); % default=1e-6
x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
        end
        eCtrd(I,II) = nanmean(jwd1); eSprd(I,II) = nanmean(jwd2);
        dCtrd(I,II) = nanmean(jwd3); dSprd(I,II) = nanmean(jwd4);
        subplot(NY,NX,(NY-I)*NX+II), hold on
        plot(f,w,'k-'), plot(stimF,r,'ro--')
        text(max(xlim),max(ylim),sprintf('Ctrd = %.2f (%.2f) kHz\nSprd = %.2f (%.2f) kHz',Ctrd,eCtrd(I,II),Sprd,eSprd(I,II)),'HorizontalALignment','right','VerticalAlignment','top')
    end
end

% fun = @(x,xdata)1+(x(1)*x(2)).*exp(-x(1)*x(2))*xdata

figure(2), clf
colormap('gray')
subplot(1,2,1), contourf(C,NStim,dCtrd), set(gca,'XScale','log','YScale','log'), title('Ctrd'), xlabel('C'), ylabel('NStim')
subplot(1,2,2), contourf(C,NStim,dSprd), set(gca,'XScale','log','YScale','log'), title('Sprd'), xlabel('C')
end
% Local functions
function nErb = lcfF2NErb(f)
nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);
end
function f = lcfNErb2F(nErb)
f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
end
function erb = lcfErb(f)
erb = 24.67*(4.37*f+1)/1000;
end
function [ECtrd,ESprd,r] = lcfVoxResp(CF,C,NOI,f)
ERB = C*lcfErb(CF);
P = 4*CF/ERB;  
g = abs(f-CF)/CF;
r = (1+P*g).*exp(-P*g)+NOI/100*randn(size(f));
ECtrd = sum(max(r,0).*f)/sum(max(r,0));
ESprd = sqrt(sum(max(r,0).*(f-ECtrd).^2)/sum(max(r,0)));
end
function [Ctrd,Sprd,f,w] = lcfW(CF,C)
ERB = C*lcfErb(CF);
P = 4*CF/ERB; 
f = CF+(-3*CF:(6*CF)/(100-1):3*CF); f = f(f>=0);
g = abs(f-CF)/CF;
w = (1+P*g).*exp(-P*g);

Ctrd = sum(w.*f)/sum(w);
Sprd = sqrt(sum(w.*(f-Ctrd).^2)/sum(w));
end


function [f,w] = lcfROEX(param,f)
C =param(1);
CF = param(2);
ERB = C*lcfErb(CF);
P = 4*CF/ERB; 
g = abs(f-CF)/CF;
w = (1+P*g).*exp(-P*g);
end

function hrf = makeHrf(TR) 
% given the TR, return the HRF shape for t = 0 ... 30s
%
% using the equation given in the lecture (simple boynton version)
tau = 2; % time constant (s) 
delta = 2; % time shift (s)
t = [0:TR:30]; % vector of time points (in steps of TR)
tshift = max(t-delta,0); % shifted, but not < 0
hrf = (tshift ./ tau) .^2 .* exp(-tshift ./tau )./(2*tau); 
 
end





