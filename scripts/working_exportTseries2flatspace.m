% Load time series (in volumetric space)
% get transformation info (from vol to flat space)
% apply transformation to time series
% save to new/correct group

% average over depth?
% create new group with one slice (average of desired depths)?

  tseries = loadTSeries(viewBase,scanNum);
  
  
          [overlayData(iOverlay).data{iScan}, voxelSize, baseCoordsMap{iScan}] = getBaseSpaceOverlay(thisView, overlayData(iOverlay).data{iScan},[],[],baseSpaceInterp);
  
%   ============ 'base2base' ============
   xform = viewGet(view,'base2base',[baseNum1],[baseNum2])
%    This will return the xform matrix that specifies the
%    transformation from the current base coordinates to the specified
%    base (baseNum1)'s coordinates. If you specify baseNum2 it will
%    calculate the xform matrix from baseNum2 to baseNum1

% view = saveNewTSeries(view,tseries,[scanParams],[hdr],[makeLink])
%
% tseries: tseries array (x,y,z,t) or a filename of a nifti file
%
% scanParams: structure with the following fields: fileName, description,
% junkFrames, nFrames. The other scanParams fields are set automatically
% from the nifti header. 

% look at concat scan - how they save new group and scan

% Open new view and set its group to the concat group name. Create the
% group if necessary.
viewConcat = newView;
concatGroupNum = viewGet(viewConcat,'groupNum',params.newGroupName);
if isempty(concatGroupNum)
  view = viewSet(view,'newgroup',params.newGroupName);
  concatGroupNum = viewGet(viewConcat,'groupNum',params.newGroupName);
end
viewConcat = viewSet(viewConcat,'currentGroup',concatGroupNum);

% then save new time series
  
  %save the output
  if params.exportToNewGroup
    baseName = [viewGet(thisView,'currentBaseName') 'Volume'];
    groupName=baseName;
    if ismember(groupName,viewGet(thisView,'groupNames'))
      fprintf('(combineTransformOverlays) Installing output overlays in group %s.\n',groupName);
      thisView = viewSet(thisView,'currentGroup',groupName);
      thisView = viewSet(thisView,'currentBase',viewGet(thisView,'baseNum',groupName));
    else
      fprintf('(combineTransformOverlays) Creating group %s.\n',groupName);
      thisView = viewSet(thisView,'newGroup',groupName);
      thisView = viewSet(thisView,'currentGroup',groupName);
      % export the flat map as an empty volume
      base = viewGet(thisView,'baseCache');
      baseNum = viewGet(thisView,'currentBase');
      if isempty(base)
        base.im = getBaseSlice(thisView,viewGet(thisView,'curslice'),viewGet(thisView,'baseSliceIndex',baseNum),viewGet(thisView,'rotate'),baseNum,viewGet(thisView,'basetype'));
      end
      baseVolume = viewGet(thisView,'baseVolume');
      hdr = baseVolume.hdr;
      hdr.bitpix = 32;   
      hdr.datatype = 16;
      hdr.is_analyze = 1;
      hdr.scl_slope = 1;
      hdr.endian = 'l';
      voxelSize(1:2) = repmat(mean(voxelSize(1:2)),1,2);
      if any(voxelSize ~= viewGet(thisView,'basevoxelsize',baseNum))
         hdr.pixdim = [0 voxelSize 0 0 0 0]';        % all pix dims must be specified here
         hdr.qform44 = diag([voxelSize 0]);
         hdr.sform44 = hdr.qform44;
      end
      % Copy file to the tseries directory
      tseriesDir = viewGet(thisView,'tseriesDir');
      scanFileName = [baseName mrGetPref('niftiFileExtension')];
      newPathStr = fullfile(tseriesDir,scanFileName);
      [bytes,hdr] = cbiWriteNifti(newPathStr,repmat(base.im,[1 1 size(outputData{1},3)]),hdr);
      % Add it
      scanParams.fileName = scanFileName;
      thisView = viewSet(thisView,'newScan',scanParams);