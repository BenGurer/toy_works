function runBatMult

load('IndData.mat','exp2'), exp2 = exp2(:,:,2:3); NParts = size(exp2,1);
PF = 1; df = [0 200 600 1800]; af = PF*2.^(df/1200); Lev = 60;
load('LSQSingle.mat','Al'), M = -1.044321/1000;

% Create bootstrap smpls;
if isempty(dir('BtspSmpls.mat'))
    smpls(1,:) = 1:NParts;
    NSmpls = 1000;
    I = 2; 
    while I<=NSmpls 
        smpls(I,:)= sort(round(1+(NParts-1)*rand(1,NParts)));
        smpls = unique(smpls,'rows'); 
        if size(smpls,1)==I, I = I+1; end
    end
    save('BtspSmpls.mat','smpls')
else
    load('BtspSmpls.mat')
    NSmpls = size(smpls,1);
end
AvgIdx = find(all(smpls==repmat(1:NParts,size(smpls,1),1),2));

if ~isempty(dir('BTSPMult.dat'))
   N = lcfN('BTSPMult.dat')
   FId = fopen('BTSPMult.dat','at');
else
   N = 0
   FId = fopen('BTSPMult.dat','wt');
   fprintf(FId,'Ctip Ctail W RMSD\n');
end

for I = N+1:NSmpls
% for I = AvgIdx:AvgIdx
    y = squeeze(nanmean(exp2(smpls(I,:),:,:)));
    if ~any(isnan(y(:)))   
        [Ctip,Ctail,W,RMSD,ad2,ad3] = batLSQMult(PF,af,Lev,Al,abs(M),y);
        fprintf(FId,'%d %g %g %g %g\n',I,Ctip,Ctail,W,RMSD);

        figure(3), clf, hold on
        plot(df,y(:,1),'ks-'), plot(df,y(:,2),'ko--')
        plot(df,ad2,'rs-'), plot(df,ad3,'ro--')
        set(gca,'YScale','log')
        text(max(xlim),max(ylim),sprintf('%d: RMSD = %g%%',I,RMSD),'HorizontalAlignment','right','VerticalAlignment','top')            
    else
        fprintf(FId,'%d %g %g %g %g\n',I,nan,nan,nan,nan);
    end              
end
fclose(FId);

% ***** lcfN *****
function N = lcfN(file)

FId = fopen(file);
fgetl(FId);
while ~feof(FId)
    [tok,rest] = strtok(fgetl(FId));
    N = sscanf(tok,'%d');
end
fclose(FId);


