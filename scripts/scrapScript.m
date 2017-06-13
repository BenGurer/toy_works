%GLM analysis
thisView = viewSet(thisView,'curGroup','ConcatenationHLsim');
[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.hrfModel = 'hrfBoxcar';
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_BoxCar';
glmParams.hrfParams.description = 'GLM Box Car -HLsim Concat';
glmParams.hrfParams.delayS =  2.5;
glmParams.hrfParams.durationS = 2.5;
glmParams.scanParams{1}.stimDurationMode = 'fromFile';
glmParams.scanParams{1}.supersamplingMode =  'Set value';
glmParams.scanParams{1}.designSupersampling = 3;
glmParams.scanParams{1}.acquisitionDelay = .75;
glmParams.numberFtests = 1;
glmParams.fTestNames{1} = 'fTest 1';
glmParams.numberContrasts = 0;
glmParams.parametricTests = 0;
glmParams.outputEstimatesAsOverlays = 1; 
[thisView, glmParams] = glmAnalysis(thisView,glmParams);

%Tonotopy analysis
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);

params.combineFunction='weightedMeanStd';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);

%% save analysis
saveAnalysis(thisView,'GLM_BoxCar')