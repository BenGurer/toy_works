function runBatSingle

% Load individual data; this function fits the Experiment-1 data; the data
% are adaptation in percent [(A-P)/A*100 = (1-P/A)*100, where A = unadapter response, P =
% adapted response as a function of the adapter-probe frequecy difference
% (df in cent) and adapter-probe SOA.
load('IndData.mat','exp1'), NParts = size(exp1,1); 
PF = 1; df = [0 200 600 1800]; af = PF*2.^(df/1200); 
soa = [125 250 500 1000]; Lev = 60;

% Create bootstrap samples; we are fitting bootstrap averages of our data.
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

% This creates a file with the model fit parameters for each bootstrap
% sample.
if ~isempty(dir('BTSPSingle.dat'))
   N = lcfN('BTSPSingle.dat')
   FId = fopen('BTSPSingle.dat','at');
else
   N = 0
   FId = fopen('BTSPSingle.dat','wt');
   fprintf(FId,'Ctip Ctail W Al M RMSD\n');
end

%Look through the samples.
for I = N+1:NSmpls
% for I = AvgIdx:AvgIdx
    y = squeeze(nanmean(exp1(smpls(I,:),:,:)));
    if ~any(isnan(y(:)))
        [M,y0] = funRegr(soa,log(y)); y0 = exp(y0); % this is the decay rate of adaptation with increasing SOA.
        % This is where the parameters of the tonotopic population model
        % are fitted. If fitting is appropriate in your case, it might be to
        % individual neuron responses or population averages.
        % Alternatively, a representative model (not fitted to any data
        % that shows the relevant features might be better. In our model, we
        % fitted the neuron frequency response functions as a two-filter response, 
        % with a tip and a broader tail. Ctip and Ctail are factor that describe the
        % widths of eth tip and tail filetrs in relation to the Glasberg and Moore 
        % (1990) normal ERB auditory filters (humans). Again, a simple roex function 
        % might be sufficient in your case.  
        [Ctip,Ctail,W,Al,RMSD,ad] = batLSQSingle(PF,af,Lev,soa,abs(M),y,y0); 
        fprintf(FId,'%d %g %g %g %g %g %g\n',I,Ctip,Ctail,W,Al,M,RMSD);

        figure(1), clf
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
        fprintf(FId,'%d %g %g %g %g %g %g\n',I,nan,nan,nan,nan,nan,nan);
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
    




    

