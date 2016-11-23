function alignFunctional2HighResT2Star(xformFile,HighResT2Star,filenames)
% We assume the s-form in this header transforms from high res t2Star to whole head space
T2StarHDR = mlrImageReadNiftiHeader(HighResT2Star);
% shift0p5 = [1 0 0 0.5; 0 1 0 0.5; 0 0 1 0; 0 0 0 1];

% xform = transformation from the fMRI to the low resolution (same as fMRI)
% high res t 2 star - This is why we have to take into account the difference in in-plane voxel size
xform = textread(xformFile);
xform = xform(1:4,1:4);
if ieNotDefined('filenames')
    % Prompt user for filename(s)
    filenames = mlrGetPathStrDialog(pwd,'Choose one or more nifti files',{'*.img;*.nii','NIFTI Files'},'on');
end
if ~iscell(filenames)
    filenames = {filenames};
end
for p = 1:length(filenames)
    hdr = mlrImageReadNiftiHeader(filenames{p});
    voxelSizeRatio = T2StarHDR.dim(2:3)./hdr.dim(2:3);
    fMRI2highResT2Star = [voxelSizeRatio(1) 0 0 0; 0 voxelSizeRatio(2) 0 0;0 0 1 0; 0 0 0 1];
    
    % compose transform matrix
    % Transform from fMRI to low res T2 star (xform), shift voxel index to centre (shift0p5), then convert from low
    % resolution to high resolution voxel size (fMRI2highResT2Star), shift
    % back voxel index, final, apply transformation from high res t2 star to whole head
    %     hdr = cbiSetNiftiSform(hdr,T2StarHDR.sform44 * inv(shift0p5) * fMRI2highResT2Star * shift0p5 * xform);
    
    % compose transform matrix
    % Transform from fMRI to low res T2 star (xform), convert from low resolution to high resolution voxel size (fMRI2highResT2Star) and apply
    % transformation from high res t2 star to whole head
    hdr = cbiSetNiftiSform(hdr,T2StarHDR.sform44* fMRI2highResT2Star * xform);
    
    hdr.sform_code = 1;
    hdr = mlrImageWriteNiftiHeader(hdr,filenames{p});
end
