function testSglShp

[~,sglRMSD] = lcfReadDat('BTSPSingle.dat');
[~,shpRMSD] = lcfReadDat('BTSPSharp.dat');

figure(1), subplot(2,1,1) 
hist(sglRMSD,25)
set(gca,'XLim',[0 10])
title('Single')
subplot(2,1,2) 
hist(shpRMSD,25)
set(gca,'XLim',[0 10])
title('Sharp')

[H,P,ks2stat] = kstest2(sglRMSD,shpRMSD)

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
    
