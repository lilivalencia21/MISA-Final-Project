% Probabilistic Atlas 
% Gülnur Ungan - Liliana Valencia
% MISA Final Project 
% MAIA2 2018

clc; clear all; close all

%Load the first image: use the header to save the resulting image
% To save the label's probabilities 
Filename_labels = fullfile('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\TrainingRegisteredLabels\IBSR_01_seg\result.nii.gz');

label = load_untouch_nii(Filename_labels);
labels = label.img;
%To save the T1 Template
Filename_intensity = fullfile('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\TrainingRegisteredImg\IBSR_01\result.1.nii.gz');
intensity = load_untouch_nii(Filename_intensity);
int = intensity.img;

% Define paths 
% filenames = dir('D:\MAIA\Spain\MISA\Project\TrainingRegisteredLabels');
% filenames2 = dir ('D:\MAIA\Spain\MISA\Project\TrainingRegisteredImg');
filenames = dir('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\TrainingRegisteredLabels');
filenames2 = dir ('D:\MAIA\Spain\MISA\ProjectRegisfirstImage\TrainingRegisteredImg');


%Initialize parameters 
numOfFiles = length(filenames);
atlas_CSF = [];
atlas_GM = [];
atlas_WM = [];
T1 = [];

% Split images according to the tissue type and concatenate the matrices 
for fileidx = 3 : numOfFiles
    filename = strcat(filenames(fileidx).folder,"\",filenames(fileidx).name,"\result.nii.gz");
    filename = char(strrep(filename,"\","/"));
    fprintf("Filename = %s\n",filename);
    img = load_untouch_nii(filename);
    img = img.img;
    fprintf("image loaded\n");
    
    filename2 = strcat(filenames2(fileidx).folder,"\",filenames2(fileidx).name,"\result.1.nii.gz");
    filename2 = char(strrep(filename2,"\","/"));
    fprintf("Filename = %s\n",filename2);
    img_intensity = load_untouch_nii(filename2);
    img_intensity = img_intensity.img;
    fprintf("image loaded\n");
    
    CSF = img == 1;
    GM = img == 2;
    WM = img == 3;
   
    atlas_CSF = cat(4,CSF,atlas_CSF);
    atlas_GM = cat(4,GM,atlas_GM);
    atlas_WM = cat(4,WM,atlas_WM);
    
    T1 = cat(4,T1,img_intensity);
end

%Compute the mean 
atlas_CSF_mean = mean(atlas_CSF,4);
% imshow(atlas_CSF_mean(:,:,130),[])
atlas_GM_mean = mean(atlas_GM,4);
% imshow(atlas_GM_mean(:,:,130),[])
atlas_WM_mean = mean(atlas_WM,4);
T1_template = mean(T1,4); 

%Rotate images 180 degrees before to save
% atlas_CSF_mean = imrotate(atlas_CSF_mean,180);
% atlas_GM_mean = imrotate(atlas_GM_mean,180);
% atlas_WM_mean = imrotate(atlas_WM_mean,180);
% T1_template = imrotate(T1_template,180);
% nbins = 1000;
% histogram(atlas_CSF_mean,nbins);

%Save the results
label.img = atlas_CSF_mean;
save_untouch_nii(label, 'atlas_CSF.nii.gz')
label.img = atlas_GM_mean;
save_untouch_nii(label, 'atlas_GM.nii.gz')
label.img = atlas_WM_mean;
save_untouch_nii(label, 'atlas_WM.nii.gz')
intensity.img = T1_template;
save_untouch_nii(intensity, 'T1_template.nii.gz')

