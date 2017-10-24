function data = script_averageAcrossDepths
a = viewGet(thisView,'Overlay','Scan1(Tone 100Hz,0)');


  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' num2str(curOverlay-1)]);
  params.combineFunction='averageDepthVol';
  [thisView,params] = combineTransformOverlays(thisView,params);