function testSglMlt

sglPars = lcfReadDat('BTSPSingle.dat');
mltPars = lcfReadDat('BTSPMult.dat');
erbs = funCalErb(sglPars(:,1),sglPars(:,2),sglPars(:,3));
erbm = funCalErb(mltPars(:,1),mltPars(:,2),mltPars(:,3));

figure(1), subplot(2,1,1) 
hist(sglPars(:,1),100)
% set(gca,'XLim',[0 10])
title('Single Ctip')
subplot(2,1,2) 
hist(mltPars(:,1),100)
% set(gca,'XLim',[0 10])
title('Mult Ctip')

figure(2), subplot(2,1,1) 
hist(sglPars(:,2),100)
% set(gca,'XLim',[0 10])
title('Single Ctail')
subplot(2,1,2) 
hist(mltPars(:,2),100)
% set(gca,'XLim',[0 10])
title('Mult Ctail')

figure(3), subplot(2,1,1) 
hist(sglPars(:,3),100)
% set(gca,'XLim',[0 10])
title('Single W')
subplot(2,1,2) 
hist(mltPars(:,3),100)
% set(gca,'XLim',[0 10])
title('Mult W')

figure(4), subplot(2,1,1) 
hist(erbs,100)
% set(gca,'XLim',[0 10])
title('Single ERB')
subplot(2,1,2) 
hist(erbm,100)
% set(gca,'XLim',[0 10])
title('Mult ERB')

[H,P,ks2stat] = kstest2(erbs,erbm)

% ***** lcfReadDat *****
function [pars,rmsd] = lcfReadDat(file)

FId = fopen(file);

line = fgetl(FId); nams = {};
while ~isempty(line)
    [tok,line] = strtok(line);
    nams = [nams {tok}];
end
N = length(nams);

while ~feof(FId)
    [tok,rest] = strtok(fgetl(FId));
    I = sscanf(tok,'%d'); jwd = sscanf(rest,'%g ',[1 N]);
    pars(I,:) = jwd(1:end-1); rmsd(I) = jwd(end);
end
fclose(FId);