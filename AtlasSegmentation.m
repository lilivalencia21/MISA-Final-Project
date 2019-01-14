% Atlas Segmentation 
% MISA Atlas Segmentation Project
% Gülnur Ungan & Liliana Valencia

clear all
close all
clc
    
tic
% Load Mask to remove background
filename ='D:\MAIA\Spain\MISA\ProjectRegisfirstImage\ValidationIntensityNormalized\IBSR_17_Norm.nii.gz';

nii = load_untouch_nii(filename);
mask = nii.img;
[rows,cols,idx]=size(mask);

% Load the atlas of each tissue type
filename_CSF = fullfile('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\ValAtlasRegistration\IBSR_17\CSF\result.nii.gz');
CSF = load_untouch_nii(filename_CSF);
img_CSF = CSF.img;
filename_GM = fullfile('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\ValAtlasRegistration\IBSR_17\GM\result.nii.gz');
GM = load_untouch_nii(filename_GM);
img_GM = GM.img;
filename_WM = fullfile('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\ValAtlasRegistration\IBSR_17\WM\result.nii.gz');
WM = load_untouch_nii(filename_WM);
img_WM = WM.img;

% Remove the background
mask = reshape(mask,[],1,1);
nonzeroIdx = find(mask);
zerosIdx = find(~mask);

img_CSF = reshape(img_CSF,[],1,1);
img_CSF = img_CSF(nonzeroIdx,1);
img_GM = reshape(img_GM,[],1,1);
img_GM = img_GM(nonzeroIdx,1);
img_WM = reshape(img_WM,[],1,1);
img_WM = img_WM(nonzeroIdx,1);
Atlas = [img_CSF img_GM img_WM];

% Compute the labels
[~,label] = max(Atlas');
label = label';

% Create an empty vector for the segmentation result
segmented_nonzero = zeros(length(nonzeroIdx),1);
segmented_zero = zeros(length(zerosIdx),1);
segmentation = vertcat(segmented_nonzero,segmented_zero);

% Assign the labels to the correct positions
for i = 1:length(label)
    t = nonzeroIdx(i);
    segmentation(t) = label(i);
end

% Load label image
filename_GT = 'D:\MAIA\Spain\MISA\ProjectRegisfirstImage\TrainingValidationTestSets\Validation_Set\IBSR_17\IBSR_17_seg.nii.gz';
nii_GT = load_untouch_nii(filename_GT);
GT = nii_GT.img;
GT = double(GT);
GT = reshape(GT,[],1,1);

% Compute dice score
similarity = dice(segmentation,GT);

% Save Image
segmentation = reshape(segmentation,[rows,cols,idx]);
nii.img = segmentation;
save_untouch_nii(nii,'segmentation_IBSR_17.nii.gz');
toc