function subjectInfo = getSubjectInfo_sHL(iSub)
    %
    %   usage: getSubjectInfo_sHL(iSub)
    %      by: Ben Gurer
    %    date: 18/10/2017
    % purpose: get subject information for simulated hearing loss
    %   input: subject number (iSub)
    %
    

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

%% output selected subjects daya
subjectInfo = struct();
subjectInfo.subjectID = subject{iSub};
subjectInfo.freeSurferName = freeSurferName{iSub};
subjectInfo.flatmapNames = flatmapName{iSub};