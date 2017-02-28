iSubj = 1;

epiDims = [128 128 24 73]; % dims of functional scans

% add check from startup

if ispc
    dataDir = 'N:/data';
elseif isunix
    dataDir = '/home/beng/data';
end
studyDir = 'CorticalMagnification';

sides = {'left','right'};
Sides = {'Left','Right'};

% Subject info
subjects{1} = '02344_034';
niftiBaseName{1} = [];
wholeheadMPRAGE{1} = 'Sparse_02344_MPRAGE_SENSE_13_1';
freeSurferName{1} = '02344_034';
% T2star{1} = '16';
% refScan{1} = '15'; % scan before t2 structural

% distCorrectionRefSparse{1} = {'17','18'};
% distCorrectionRefCont{2} = {'20','21'};
% % scanList = [13,15,19,22];
% distCorrectionRef = {'20','21'};

cd([dataDir '/scanner/' subjects{iSubj}])

    %% check for filenames below 10 and add a zero before the number - make sure linking of stim files works later
    % str =
    % [token, remain] = strtok(str, ...)
scanFiles = dir;
for id = 1:length(scanFiles)
    str = scanFiles(id).name;
    strParts = strsplit(str,'_');
    if length(strParts) > 1
    if strcmp([char(strParts(2)) '_' char(strParts(3))],subjects{iSubj})
        checkScanNum = char(strParts(4));
        numStrParts = length(strParts);
        if 10>str2num(checkScanNum) && numel(checkScanNum)==1
            newName = [char(strParts(1)) '_' char(strParts(2)) '_' char(strParts(3)) '_' ['0' checkScanNum] '_' char(strParts(5))];
            movefile(scanFiles(id).name,newName);
        end
    end
    end
end
    
 
!ptoa -f -q -nii *.PAR


if ~isempty(wholeheadMPRAGE{iSubj})
    % Move whole head PSIR
    mkdir(fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    movefile(fullfile(dataDir,'scanner',subjects{iSubj},[wholeheadMPRAGE{iSubj} '*.nii']),fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    cd(fullfile(dataDir,'Anatomy/originals/',freeSurferName{iSubj}));        
    %RUN recon-all DIRECTLY IN TERMINAL
%     fprintf(['recon-all_highres ' freeSurferName{iSubj} ' ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' freeSurferName{iSubj} '_PSIR_pos_-.7_thr\n'])
    %     recon-all -all -s $SUBJECT -hires -i $IMAGE -expert $EXPERT_FILE
    system(sprintf('fslmaths %s_PSIR_pos_-.7_thr -mul 200 %s_PSIR_pos_-.7_thr_200',freeSurferName{iSubj},freeSurferName{iSubj}));
    % fprintf(['recon-all -all -s ' freeSurferName{iSubj} ' -hires -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' freeSurferName{iSubj} '_PSIR_pos_-.7_thr.nii' ' -expert ' fullfile(dataDir,'Anatomy','freesurfer','high_res_options.txt\n')])
    fprintf(['recon-all -all -s ' freeSurferName{iSubj} ' -hires -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' freeSurferName{iSubj} '_PSIR_pos_-.7_thr_200.nii' ' -expert ' fullfile(dataDir,'Anatomy','freesurfer','high_res_options.txt\n')])
    fprintf(['recon-all -subjid ' freeSurferName{iSubj} ' -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' wholeheadMPRAGE{iSubj} '.nii' ' -all']);
end