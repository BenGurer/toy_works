thisView = getMLRView;
for iSide=1:2
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'baseNum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat']));
  rotation(iSide) = viewGet(thisView,'rotate');
end
params.path = fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj},'surfRelax');
for iFlat =1:5
  for iSide=1:2
    params.anatFileName = [freeSurferName{iSubj} '_mprage_pp.nii'];
    params.flatRes=3;
    params.threshold = 1;
    params.baseValues = 'curvature';
    switch(iFlat)
      case 1
        whichSurf = ['_invFNIRT_' subjects{iSubj}];
        params.baseValues = 'curvature';
      case 2
        whichSurf = ['_invFNIRT_' subjects{iSubj}];
        params.flatRes=1;
      case 3
        whichSurf = '';
        params.threshold = 0;
      case 4
        whichSurf = '';
        params.baseValues = 'thickness';
        params.threshold = 0;
      case 5
        %load subject flat sampling MNI single subject space
        params.anatFileName = fullfile(dataDir,'Anatomy/freesurfer/subjects/Colin27/surfRelax','Colin27_mprage_pp.nii');
        whichSurf = '_Colin27';
    end
    params.flatFileName = [freeSurferName{iSubj} '_' sides{iSide} '_Flat.off'];
    params.outerCoordsFileName = [freeSurferName{iSubj} '_' sides{iSide} '_GM' whichSurf '.off'];
    params.innerCoordsFileName = [freeSurferName{iSubj} '_' sides{iSide} '_WM' whichSurf '.off'];
    params.curvFileName = [freeSurferName{iSubj} '_' sides{iSide} '_Curv.vff'];
    base = importFlatOFF(params);
    switch(iFlat)
      case 1
        base.name = [freeSurferName{iSubj} '_' sides{iSide} '_Flat' whichSurf];
      case 2
        base.name = [freeSurferName{iSubj} '_' sides{iSide} '_Flat' whichSurf '_lowres'];
      case 3
        base.name = [freeSurferName{iSubj} '_' sides{iSide} '_Flat' whichSurf '_curvature'];
      case 4
        base.name = [freeSurferName{iSubj} '_' sides{iSide} '_Flat' whichSurf '_thickness'];
      case 5
        base.name = [freeSurferName{iSubj} '_' sides{iSide} '_Flat' whichSurf];
    end
    thisView = viewSet(thisView, 'newbase', base);
    thisView = viewSet(thisView,'rotate',rotation(iSide));
  end
end