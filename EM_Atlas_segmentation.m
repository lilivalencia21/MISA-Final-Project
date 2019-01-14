%% EM with Atlas initialization
% LAB 3 MISA
% Gülnur Ungan - Liliana Valencia
clc; close all; clear all
tic
% Load intensity image 
filename = 'D:\MAIA\Spain\MISA\ProjectRegisfirstImage\ValidationIntensityNormalized\IBSR_17_Norm.nii.gz';
nii_intensity = load_untouch_nii(filename);
intensity = nii_intensity.img;
[rows,cols,imgIdx]=size(intensity);


%Load the result of the atlas segmentation 
filename = 'D:\MAIA\Spain\MISA\ProjectRegisfirstImage\AtlasSegmentationResults\segmentation_IBSR_17.nii.gz';
nii_atlasSegmentation = load_untouch_nii(filename);
Atlas_segm = nii_atlasSegmentation.img;
% Atlas_segm = double(Atlas_segm);

% Remove the background
mask = double(intensity);
mask = reshape(mask,[],1,1);
% nonzeroIdx = find(mask);
% zerosIdx = find(~mask);
%Remove skull and background
k=1;
j=1;
for i=1:length(mask)
    if (mask(i)<= 0)
        %Assign a value that's not in the matrix for zero position
        zerosIdx(k,1) = (i);   
        k=k+1;
    else
        %Keep the indexes to reconstruct the result
        nonzeroIdx(j,1) = (i);     
        j=j+1;
    end
end

intensity = reshape(intensity,[],1,1);
X = intensity(nonzeroIdx,1);
X = double(X);

img_segm = reshape(Atlas_segm,[],1,1);
img_segm = img_segm(nonzeroIdx,1);
img_segm = double(img_segm);

%% Initalize parameters
%Covariance
% sigma=zeros(2,2,3);
sigma(1,1) = cov(X(img_segm==1));      
sigma(2,1) = cov(X(img_segm==2));
sigma(3,1)= cov(X(img_segm==3));
%mean
miu(1,1) = mean(X(img_segm==1));
miu(2,1) = mean(X(img_segm==2));
miu(3,1) = mean(X(img_segm==3));
%Alpha
alpha = [1/3 1/3 1/3];  

K = 3;      %Number of clusters or tissue types
N = length(X);
W = zeros(N,K);

%Computation of Weights 
gm = gmdistribution(miu(1,1),sigma(1,1));
W(:,1)=pdf(gm,X)*alpha(1);
gm = gmdistribution(miu(2,1),sigma(2,1));
W(:,2)=pdf(gm,X)*alpha(2);
gm = gmdistribution(miu(3,1),sigma(3,1));
W(:,3)=pdf(gm,X)*alpha(3);

for i=1:N
    W(i,:)= W(i,:)/sum(W(i,:)); 
end

[v,l]= max(W');
weigths = v';
label = l';

iter=7;
difference=0;
%Initialize vector for likelihood
% likelihood = zeros((iter+1),1);
% likelihood(1,1)=0;

for num_iter=1:iter
    %M-Step
    N_k = sum(W);
    alpha_new = N_k./N;

    %Compute the new mean
    
    secondterm = sum(X.*W(:,1));
    miu_new(1,:)=(1/N_k(1)).*secondterm;
    secondterm = sum(X.*W(:,2));
    miu_new(2,:)=(1/N_k(2)).*secondterm;
    secondterm = sum(X.*W(:,3));
    miu_new(3,:)=(1/N_k(3)).*secondterm;
    

    %Compute the new covariance
    sigma1 = ((W(:,1)'.*(X-miu_new(1,1))')*(X- miu_new(1,1)));
    sigma2 = ((W(:,2)'.*(X-miu_new(2,1))')*(X- miu_new(2,1)));
    sigma3 = ((W(:,3)'.*(X-miu_new(3,1))')*(X- miu_new(3,1)));
    sigma_new(1,1) = (1/N_k(1)).*sigma1;
    sigma_new(2,1) = (1/N_k(2)).*sigma2;
    sigma_new(3,1) = (1/N_k(3)).*sigma3;
      
    %E-Step
    %Compute the new weights
    gm = gmdistribution(miu_new(1,1),sigma_new(1,1));
    W(:,1)=pdf(gm,X)*alpha_new(1);
    gm = gmdistribution(miu_new(2,1),sigma_new(2,1));
    W(:,2)=pdf(gm,X)*alpha_new(2);
    gm = gmdistribution(miu_new(3,1),sigma_new(3,1));
    W(:,3)=pdf(gm,X)*alpha_new(3);
    
     
    for i=1:N
        W(i,:)= W(i,:)/sum(W(i,:)); 
    end

    %Check and assign weights and labels
    [v,l]= max(W');
    v =v';
    l = l';
      for i=1:length(W)   
        if v(i) > weigths(i)   %Compare the weight with the previous one 
            weigths(i) = v(i);
            label(i) = l(i);
        end
      end 


   
end


%% Load label image
filename_GT = 'D:\MAIA\Spain\MISA\ProjectRegisfirstImage\TrainingValidationTestSets\Validation_Set\IBSR_17\IBSR_17_seg.nii.gz';
nii_GT = load_untouch_nii(filename_GT);
GT = nii_GT.img;
GT = double(GT);
GT = reshape(GT,[],1,1);

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

img_CSF = reshape(img_CSF,[],1,1);
img_CSF = img_CSF(nonzeroIdx,1);
img_GM = reshape(img_GM,[],1,1);
img_GM = img_GM(nonzeroIdx,1);
img_WM = reshape(img_WM,[],1,1);
img_WM = img_WM(nonzeroIdx,1);

Atlas = [img_CSF img_GM img_WM];

new_W = W.*Atlas;
[~,final_label]= max(new_W');
final_label = final_label';

segmentation = zeros(length(GT),1);

% Assign the labels to the correct positions
for i = 1:length(label)
    t = nonzeroIdx(i);
    segmentation(t) = final_label(i);
end


% Compute dice score
similarity = dice(segmentation,GT);

% Save Image
segmentation = reshape(segmentation,[rows,cols,imgIdx]);
nii_GT.img = segmentation;
save_untouch_nii(nii_GT,'segmentation_EM_Atlas_IBSR_17.nii.gz')
toc