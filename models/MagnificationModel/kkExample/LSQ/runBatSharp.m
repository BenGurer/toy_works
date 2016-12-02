function runBatSharp

load('IndData.mat','exp1'), NParts = size(exp1,1);
PF = 1; df = [0 200 600 1800]; af = PF*2.^(df/1200); 
soa = [125 250 500 1000]; Lev = 60;

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

if ~isempty(dir('BTSPSharp.dat'))
   N = lcfN('BTSPSharp.dat')
   FId = fopen('BTSPSharp.dat','at');
else
   N = 0
   FId = fopen('BTSPSharp.dat','wt');
   fprintf(FId,'Ctip Ctail W Al s(1) s(2) s(3) s(4) RMSD\n');
end

for I = N+1:NSmpls
% for I = AvgIdx:AvgIdx
    y = squeeze(nanmean(exp1(smpls(I,:),:,:)));
    if ~any(isnan(y(:)))
        [~,y0] = funRegr(soa,log(y)); y0 = exp(y0);
        [Ctip,Ctail,W,Al,s,RMSD,ad] = batLSQSharp(PF,af,Lev,y,y0);
        fprintf(FId,'%d %g %g %g %g %g %g %g %g %g\n',I,Ctip,Ctail,W,Al,s(1),s(2),s(3),s(4),RMSD);

        figure(2), clf
        subplot(1,2,1), hold on
        for II = 1:length(soa)
            plot(df,y(:,II),'ks-'), axis tight
            plot(df,ad(:,II),'ro--'), axis tight
        end
        set(gca,'YScale','log')
        subplot(1,2,2), hold on
        for II = 1:length(af)
            plot(soa,y(II,:),'ks-'), axis tight
            plot(soa,ad(II,:),'ro--'), axis tight
        end
        set(gca,'YScale','log')
        text(max(xlim),max(ylim),sprintf('%d: RMSD = %g%%',I,RMSD),'HorizontalAlignment','right','VerticalAlignment','top')            
    else
        fprintf(FId,'%d %g %g %g %g %g %g %g %g %g\n',I,nan,nan,nan,nan,nan,nan,nan,nan,nan);
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

    

