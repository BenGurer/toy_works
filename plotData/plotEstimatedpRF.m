function plotEstimatedpRF
% (filenames)
% filenames={'SparseROIbetas.mat','ContROIbetas.mat'};
% filenames={'Concatenation Sparse_pRFEst_9bins.mat','Concatenation Cont_pRFEst_9bins.mat'};
dataDir = 'N:\data\CorticalMagnification\03644_012\';
filenames={'Concatenation Sparse_PAC_pRFEst_9bins.mat','Concatenation Cont_PAC_pRFEst_9bins.mat'};
dataNames = {'Sparse','Cont'};
ComprfModel= {};
CompnVoxelsPerBin= [];


figure;


for i = 1:length(filenames)
    
    load([dataDir filenames{i}])
    ComprfModel{i} = rfModel;
    CompnVoxelsPerBin(:,i) = nVoxelsPerBin;
    CompscaleMean(:,i) = scaleMean;
end
%% make legend lables
ROIbetaslabel = cell(size(nVoxelsPerBin,1),1);
for i = 1:length(nVoxelsPerBin)
    % ROIbetaslabel{i} = ['beta = ' num2str(i) ' AvSte = ' num2str(ROISteAv(i))];
    
    ROIbetaslabel{i} = ['beta = ' num2str(i-1)];
end
for i = 1:length(filenames)
    % subplot(length(filenames),2,i)
    subplot(length(filenames),4,4*(i-1)+1)
    bar(0:length(nVoxelsPerBin)-1,CompnVoxelsPerBin(:,i));
    title([dataNames{i} ' number of voxels'])
    axis tight
    ylim ([0 max(max(CompnVoxelsPerBin))]);
    % ylim ([-0.5 3.5])
    % ylim ([-2 6])
    % % xlim ([0 length(normalisedBetas(:,i))+1])
    % hold on
    % plot([0 length(betas(:,i))+1],[0 0],'--k')
    % subplot(length(filenames),2,i+length(filenames))
    subplot(length(filenames),4,4*(i-1)+2)
    bar(0:length(nVoxelsPerBin)-1,CompscaleMean(:,i));
    axis tight
    ylim ([min(min(CompscaleMean)) max(max(CompscaleMean))]);
    title('Scaling')
    
    % subplot(length(filenames),4,i+1+length(filenames))
    subplot(length(filenames),4,4*(i-1)+3)
    PLOTprfModel = ComprfModel{i};
    for n = 1:length(PLOTprfModel)
        plot(PLOTprfModel{n},'LineWidth',2)
        set(gca,'ColorOrder',jet(length(nVoxelsPerBin)))
        hold on
        axis tight
    end
    
     title('Estimated pRF')
    
    % subplot(length(filenames),2,i+2+length(filenames))
    subplot(length(filenames),4,4*(i-1)+4)
    for n = 1:length(PLOTprfModel)
        scalepRFModel{n} = PLOTprfModel{n}.*CompscaleMean(n,i)
        plot(scalepRFModel{n},'LineWidth',2)
        set(gca,'ColorOrder',jet(length(nVoxelsPerBin)))
        hold on
        axis tight
    end
    title('Scaled Estimated pRF')
    legend(ROIbetaslabel,'Location','southeastoutside')
    % title(['mean Ste = ' num2str(mean(normalisedSte(:,i)))])
    % ylim ([-1 2.5])
    % ylim ([-1.5 3])
    % xlim ([0 length(normalisedBetas(:,i))+1])
    
end


%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,fitParams)
    
    p.x = params(1);
    %         p.y = params(2);
    p.y = 1;
    p.std = params(3);
    % use a fixed single gaussian
    p.canonical.type = 'gamma';
    p.canonical.lengthInSeconds = 25;
    p.canonical.timelag = fitParams.timelag;
    p.canonical.tau = fitParams.tau;
    p.canonical.exponent = fitParams.exponent;
    p.canonical.offset = 0;
    p.canonical.diffOfGamma = fitParams.diffOfGamma;
    p.canonical.amplitudeRatio = fitParams.amplitudeRatio;
    p.canonical.timelag2 = fitParams.timelag2;
    p.canonical.tau2 = fitParams.tau2;
    p.canonical.exponent2 = fitParams.exponent2;
    p.canonical.offset2 = 0;
    if fitParams.fitHDR
        p.canonical.exponent = params(4);
        p.canonical.timelag = params(5);
        p.canonical.tau = params(6);
        p.canonical.diffOfGamma = fitParams.diffOfGamma;
        if fitParams.diffOfGamma
            p.canonical.exponent2 = params(7);
            p.canonical.amplitudeRatio = params(8);
            p.canonical.timelag2 = params(9);
            p.canonical.tau2 = params(10);
            p.canonical.offset2 = 0;            
        end
    end

%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
  fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

% exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);

%%%%%%%%%%%%%%%%%%%%
%%   getRFModel   %%
%%%%%%%%%%%%%%%%%%%%
function rfModel = getRFModel(params,fitParams)

rfModel = [];

% convert stimulus spacing to voxel magnifcation domain
if any(strcmp(fitParams.voxelScale,{'lin'}))
     x = fitParams.stimX;
    mu = params.x;
    sigma = params.std;
elseif any(strcmp(fitParams.voxelScale,{'log'}))
    x = log10(fitParams.stimX);
    mu = log10(params.x);
    sigma = params.std;
elseif any(strcmp(fitParams.voxelScale,{'erb'}))
    x = funNErb(fitParams.stimX);
    mu = funNErb(params.x);
    sigma = params.std;
else
  disp(sprintf('(pRFFit:getRFModel) Unknown voxelScale: %s',fitParams.voxelScale));
end


% now gernerate the rfModel
if any(strcmp(fitParams.rfType,{'gaussian'}))
    rfModel = makeRFGaussian(params,fitParams,x,mu,sigma);
elseif any(strcmp(fitParams.rfType,{'ROEX'}))
    rfModel = makeRFROEX(params,fitParams,x,mu,sigma);
else
  disp(sprintf('(pRFFit:getRFModel) Unknown rfType: %s',fitParams.rfType));
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeRFGaussian   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function rfModel = makeRFGaussian(params,fitParams,x,mu,sigma)

% compute rf
% rfModel = exp(-(((fitParams.stimX-params.x).^2)/(2*(params.std^2))+((fitParams.stimY-params.y).^2)/(2*(params.std^2))));

rfModel = exp(-(((x-mu).^2)/(2*(sigma^2))+((fitParams.stimY-params.y).^2)/(2*(sigma^2))));
