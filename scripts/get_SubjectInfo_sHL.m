function subjectInfo = getSubjectInfo_sHL(iSub)
%
%   usage: getSubjectInfo_sHL(iSub)
%      by: Ben Gurer
%    date: 18/10/2017
% purpose: get subject information for simulated hearing loss
%   input: subject number (iSub)
%

% pre-allocate
totalSubjects = 5;
wholeheadMPRAGE = cell(1,totalSubjects);

% define subject info
subject{1} = '02344_034';
niftiBaseName{1} = 'HL_02344_034_';
wholeheadMPRAGE{1} = '13';
freeSurferName{1} = '02344_034';
T2star{1} = '12';
refScan{1} = '11'; % scan before t2 structural
% flatmapName{1} = {'80_131_81_Rad60', '181_127_85_Rad60'};

flatmapName{1} = {'x02344_034_left_WM_Flat_80_131_81_Rad60', 'x02344_034_right_WM_Flat_181_127_85_Rad60'};
apScan{1} = 5;
paScan{1} = 6;
epiScans{1} = {'07', '09', '10', '11'};

subject{2} = '12013_002';
niftiBaseName{2} = 'HL_12013_002_';
wholeheadMPRAGE{2} = '11';
freeSurferName{2} = '12013_002';
T2star{2} = '8';
refScan{2} = '05'; % scan before t2 structural
flatmapName{2} = {' ', ' '};
apScan{2} = 5;
paScan{2} = 6;
epiScans{2} = [03, 04, 08, 09];
flatmapName{2} = {'x12013_002_left_WM_Flat_84_99_86_Rad75','x12013_002_right_WM_Flat_175_101_79_Rad75'};

subject{3} = '12023_002';
niftiBaseName{3} = 'HL_12023_002_';
wholeheadMPRAGE{3} = '10';
freeSurferName{3} = '12023_002';
T2star{3} = '7';
refScan{3} = '04'; % scan before t2 structural
distCorrectionRef{3} = {'5','6'};
flatmapName{3} = {'x12023_002_left_WM_Flat', 'x12023_002_right_WM_Flat'};

subject{4} = '11108_007';
niftiBaseName{4} = 'HL_11108_007_';
wholeheadMPRAGE{4} = '10';
freeSurferName{4} = '11108_007';
T2star{4} = '7';
refScan{4} = '04'; % scan before t2 structural
flatmapName{4} = {'x11108_007_left_Flat', 'x11108_007_right_Flat'};

subject{5} = '13016_001';
niftiBaseName{5} = 'HL_13016_001_';
wholeheadMPRAGE{5} = 'MPRAGE_2';
freeSurferName{5} = '13016_001';
T2star{5} = 'High_res_t2__9';
fMRIScans{5} = {'4', '5', '10', '11'};
refScan{5} = '5'; % scan before t2 structural
nScans(5) = 4;
conditionOrder{5} = {[1,3],[2,4]};
flatmapName{5} = {'x13016_001_left_WM_Flat_72_101_82_Rad55', 'x13016_001_right_WM_Flat_172_112_96_Rad55'};

%% output selected subjects daya
subjectInfo = struct();
subjectInfo.subjectID = subject{iSub};
subjectInfo.niftiBaseName = niftiBaseName{iSub};
subjectInfo.wholeheadMPRAGE = wholeheadMPRAGE{iSub};
subjectInfo.fMRIScans = fMRIScans{iSub};
subjectInfo.nScans = nScans(iSub);
subjectInfo.refScan = refScan{iSub};
subjectInfo.T2star = T2star{iSub};
subjectInfo.conditionOrder = conditionOrder{iSub};
subjectInfo.freeSurferName = freeSurferName{iSub};
subjectInfo.flatmapNames = flatmapName{iSub};