
function thisView = setupAfterLocationChange(thisView,Info,subjectInfo)
    %
    %   usage: setupAfterLocationChange(thisView,Info,subjectInfo)
    %      by: Ben Gurer
    %    date: 05/22/2018
    % purpose: import flatmaps
    %   input: mrView, Study information and subject information
    %  output: mrView with flatmap path names that match current directory
    %
    
%import flat maps
newPath = fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName,'/surfRelax');

baseCoordMap = viewGet(v,'baseCoordMap',baseNum);
baseCoordMapFields = fields(baseCoordMap);
% ~isempty(strfind(baseCoordMapFields{fieldNum},'path'))
baseCoordMap.path = [newPath baseCoordMap.name];


viewSet(v,'baseCoordMap',baseCoordMap,baseNum);

end
