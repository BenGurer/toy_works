function startup(whichMrTools)


if ~exist('whichMrTools','var')
  whichMrTools ='mrTools';
end

%set the debugger to stop if error
dbstop if error

%------------ Paths ----------------------------------%
if ispc
    mountPoint = 'N:';
    if exist(mountPoint,'dir')
    disp(['Using path from Network Drive (' mountPoint ').']);
    %disp('Remember to use ''rehash path'' if you modify files outside this session of Matlab');
    %However, ther might be a problem in that Matlab will not update the path if it's changed outside from it (e.g. from the Mac OS)
    %Let's see ...
    else
        keyboard
        %problem with the network drive...?
    end
elseif isunix
    mountPoint = '/home/beng';
else
    keyboard
end
% StartDirectory = [mountPoint '/matlab/'];
StartDirectory = mountPoint;
cd (StartDirectory);
mrToolsDirectory = [mountPoint '/matlab/' whichMrTools];
disp('adding paths...')
% image processing toolbox
if ~license('test','Image_Toolbox')
  addpath(genpath([mountPoint '/matlab/Image Processing Toolbox']));
end
% Nottingham functions 
addpath(genpath([mountPoint '/matlab/nottingham']));
disp('added Nottingham path...')
% TDT 
addpath(genpath([mountPoint '/matlab/tdtMRI']));
disp('added tdtMRI path...')
% EEG lab
% addpath(genpath([mountPoint '/matlab/eeglab12_0_2_0b']));
% bgfunctions
% addpath(genpath([mountPoint '/matlab/bg-Fun']));
% disp('added bgFun path...')
% addpath(genpath([mountPoint '/matlab/auditory-functions']));
% disp('added auditoryFun path...')
% ToyWorks functions
addpath(genpath([mountPoint '/matlab/toy_works']));
disp('added toy_works path...')
% Auditory pRF
addpath(genpath([mountPoint '/matlab/auditory-pRF']));
disp('added auditory-pRF path...')
% kkfunctions
% addpath(genpath([mountPoint '/matlab/kkfunctions']));
% disp('added kkfunctions path...')
% mrTools directory
addpath(genpath(mrToolsDirectory));
disp('added mrTools path...')
%fMRIqa
addpath(genpath([mountPoint '/matlab/qa']));
disp('added fMRIqa path...')
% %shadow functions
% addpath(genpath([mountPoint '/matlab/shadow']));
% disp('added shadow functions path...')
% if ispc
%   addpath(genpath([mountPoint '/matlab/shadowWin']));
% end
disp('adding paths...Done!')
disp(['You are using ' mrToolsDirectory]);
% addpath(genpath([mountPoint '/matlab']));

%------------mrTools--------------------
% this has to be after adding the path
setpref('mrLoadRet','mrDefaultsFilename',[mountPoint '/matlab/.mrDefaults']);
mrSetPref('volumeDirectory',[mountPoint '/data/Anatomy/freesurfer/subjects']);
mrSetPref('magnet',{'Nottingham 7T' 'Nottingham 3T'});
mrSetPref('coil',{'Sense Head 32 Channels' 'Sense Head 16 Channels' 'Sense 8 Channels' 'Flex S'});
mrSetPref('pulseSequence',{'2D Gradient Echo' '3D Gradient Echo' 'Spin Echo' '3D Flash'});
mrSetPref('defaultPrecision','single');
mrSetPref('graphWindow','Make new');
mrSetPref('niftiFileExtension','.nii');
mrSetPref('checkParamsConsistency','no');
mrSetPref('verbose','no');
mrSetPref('pluginPaths',[mountPoint '/matlab/' whichMrTools '/mrLoadRet/PluginAlt/Nottingham,' ...
                        mountPoint '/matlab/somato-pRF,'...
                        mountPoint '/matlab/auditory-pRF,'...
                        ]); 
mrSetPref('interrogatorPaths',[mountPoint '/matlab/nottingham/julien/mrToolsInterrogators']); 
mrSetPref('interrogatorPaths',[mountPoint '/matlab/toy_works/mrToolsInterrogators']); 

clear mountPoint;

if ispc

elseif isunix
  mrSetPref('fslPath','/usr/local/bin/');
  %mrSetPref('fslPath','/usr/local/fsl_old/bin/'); %I use the old version because fslmaths -tfce has a problem in the new version
  
  %set opengl to software, but only if the current function is called at startup (i.e. by martlabrc)
  stack = dbstack;
  if length(stack)>1 && strcmp(stack(2).name,'matlabrc') && verLessThan('matlab','8.3')
    opengl software
  end
  
% Set environment variables
  %------------ FreeSurfer -----------------------------%
  freeSurferHome = getenv('FREESURFER_HOME');
  if exist(freeSurferHome,'dir') ~= 7
     freeSurferHome = '/usr/local/freesurfer';
     setenv('FREESURFER_HOME',freeSurferHome);
  end
  path(path,sprintf('%s/matlab',freeSurferHome));

  setenv('MINC_BIN_DIR',[freeSurferHome '/mni/bin']);   
  setenv('MINC_LIB_DIR',[freeSurferHome '/mni/lib']);
  setenv('MNI_INSTALL_DIR',[freeSurferHome '/mni']);
  setenv('PATH', [freeSurferHome '/bin:' freeSurferHome '/mni/bin:' freeSurferHome '/mni/lib:' freeSurferHome '/mni:' getenv('PATH') ]);

  system('source /usr/local/freesurfer/SetUpFreeSurfer.csh >/dev/null'); 
  disp('done freesurfer stuff...')

  %------------ FreeSurfer FAST ------------------------%
  fsFastHome = getenv('FSFAST_HOME');
  if exist(fsFastHome,'dir') ~= 7
      fsFastHome = [freeSurferHome '/fsfast'];
  end
  path(path,sprintf('%s/toolbox',fsFastHome));
  disp('done freesurfer FAST stuff...')

  clear freeSurferHome fsFastHome;

% %   %------------ psignifit ------------------------%
% %   setenv('PATH', ['/usr/local/psignifit:' getenv('PATH') ]);
% %   
% %   %------------ AFNI ------------------
% %   setenv('PATH', ['/usr/local/afni:' getenv('PATH') ]);
% % 

  %------------ FSL -------------------------
  %reset FSL output type to NIFTI_PAIR, so our tools work w/ default file format.
  setenv('FSLOUTPUTTYPE','NIFTI');  
  fslHome = getenv('FSLDIR');
  if exist(fslHome,'dir') ~= 7
     fslHome = '/usr/local/fsl';
     setenv('FSLDIR',fslHome);
  end
  system(['source ' fslHome '/etc/fslconf/fsl.csh']);
  setenv('PATH', [fslHome '/bin:' getenv('PATH')])
  disp('done FSL stuff...')

  clear fslHome
  %-----------------------------------------------------%

end

