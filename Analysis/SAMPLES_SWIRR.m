cd('First attempt worktree\Analysis')
clear all 
close all 

load('PEN_SWIR.mat')

mean_image = mean(im_reflectance, 3);
figure, imshow(mean_image, [])
title('Mean Image Across All Wavelength Bands for Mix')

%% Creating mask through the segmentation of the samples

% standard deviation image
im_uf = unfold(im_reflectance);
im_uf_sd = std(im_uf'); 
im_sd = reshape(im_uf_sd, size(im_reflectance, 1), size(im_reflectance, 2));

%Finding the threshold and making the mask 
figure
subplot (1,2,1)
imshow(im_sd,[]);

subplot(1,2,2)
histogram(im_sd,100);
title('Pixel Intensities of PC1')

% Applying the mask to the mean image
mask_mix = im_sd>2.3;
figure, imshow(mask_mix.*mean_image,[])
title('Mean Image with Mask Applied')

%Filling in the holes in the mask 
figure,
imshow(mask_mix,[]);
mask2 = imfill(mask_mix,"holes");
figure,
imshow(mask2,[]);

%Labelling the samples
L= bwlabel(mask2);
figure, 
imshow(L,[]);
colormap('jet')

% % Displaying the 10nth label 
% figure, 
% imshow(L==10,[]);
% colormap('jet')

%% Commiting to alternative approach 
x = zeros(1,25);  %for storing x coordinates, or the column indices of the centre of each sample
y = zeros(1,25);  %for storing y coordinates, or the row indices of the centre of each sample

for i = 1:25
    fprintf('Select the central region of sample %d\n',i);
    [ xi, yi, ~] = impixel(im_reflectance(:,:,10),[]);
    x(i) = round(xi(1));  % Store the x coordinate (column index)
    y(i) = round(yi(1));  % Store the y coordinate (row index)
end

L_new = zeros(size(L));  % Initialize a new label matrix with the same size as L

for i = 1:25
    L_old = L(y(i), x(i));  % Get the old label at the selected coordinates
    L_new(L == L_old) = i;  % Replace all pixels with the old label with the new label i
end
unique_labels = unique(L_new);
disp('Unique labels in L_new:');
disp(unique_labels);

figure,
imshow(L_new > 0, []); % Display the new label matrix as a binary image
colormap('jet');
title('New Label Matrix with 25 Samples Labeled')

%% Creating the table of mean reflectance per sample
nSamples = max(L_new(:));
nBands   = size(im_reflectance, 3);

im_reshaped = reshape(im_reflectance, [], nBands);
meanSpectra = nan(nSamples, nBands);

for i = 1:nSamples
    samplePixels = im_reshaped(L_new(:) == i, :);
    samplePixels = samplePixels(any(samplePixels, 2), :);
    if ~isempty(samplePixels)
        meanSpectra(i, :) = mean(samplePixels, 1);
    end
end

sampleNames = "PEN_SWIR_Sample_" + compose("%02d", (1:nSamples)');

% One cell per sample, each cell contains a 1 x nBands mean spectrum
meanSpectrumCell = mat2cell(meanSpectra, ones(nSamples,1), nBands);

spectraSummaryTable = table( ...
    cellstr(sampleNames), ...
    meanSpectrumCell, ...
    repmat(nBands, nSamples, 1), ...
    'VariableNames', {'Sample_Name', 'Mean_Spectrum', 'Number_Of_Wavelengths'} ...
);

% Save into a separate output folder
outputDir = fullfile(pwd, 'analysis_outputs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fileOutSWIR = fullfile(outputDir, 'SWIR_Summary.mat');

% Current run table
newSWIR = table( ...
    cellstr(sampleNames), ...
    meanSpectrumCell, ...
    repmat(nBands, nSamples, 1), ...
    'VariableNames', {'Sample_Name', 'Mean_Spectrum', 'Number_Of_Wavelengths'} ...
);

% Append rows if file exists, else create
if isfile(fileOutSWIR)
    S = load(fileOutSWIR, 'SWIR_Summary');
    SWIR_Summary = S.SWIR_Summary;

    % Optional: skip duplicates by sample name
    keep = ~ismember(newSWIR.Sample_Name, SWIR_Summary.Sample_Name);
    newSWIR = newSWIR(keep, :);

    SWIR_Summary = [SWIR_Summary; newSWIR];   % vertical concatenation
else
    SWIR_Summary = newSWIR;
end

save(fileOutSWIR, 'SWIR_Summary');





















%% Merging the tables
outputDir = fullfile(pwd, 'analysis_outputs');

vnirFile = fullfile(outputDir, 'mean_sample_summary.mat'); % or change to your VNIR filename
swirFile = fullfile(outputDir, 'SWIR_Summary.mat');

if ~isfile(vnirFile)
error('VNIR file not found: %s', vnirFile);
end
if ~isfile(swirFile)
error('SWIR file not found: %s', swirFile);
end

S1 = load(vnirFile);
S2 = load(swirFile);

%% Load both tables
S1 = load(vnirFile);
S2 = load(swirFile);

T_vnir = S1.mean_sample_summary;
T_swir = S2.SWIR_Summary;

%% Convert sample names to string
T_vnir.Sample_Name = string(T_vnir.Sample_Name);
T_swir.Sample_Name = string(T_swir.Sample_Name);

%% Create matching sample ID
% Example:
% BEM_VNIR_Sample_01  ->  BEM_Sample_01
% BEM_SWIR_Sample_01  ->  BEM_Sample_01

T_vnir.Sample_ID = erase(T_vnir.Sample_Name, "_VNIR");
T_swir.Sample_ID = erase(T_swir.Sample_Name, "_SWIR");

%% Rename variables before joining
T_vnir = renamevars(T_vnir, ...
    ["Sample_Name", "Mean_Spectrum", "Number_Of_Wavelengths"], ...
    ["VNIR_Name", "VNIR_Mean_Spectrum", "VNIR_Number_Of_Wavelengths"]);

T_swir = renamevars(T_swir, ...
    ["Sample_Name", "Mean_Spectrum", "Number_Of_Wavelengths"], ...
    ["SWIR_Name", "SWIR_Mean_Spectrum", "SWIR_Number_Of_Wavelengths"]);

%% Move Sample_ID to the front
T_vnir = movevars(T_vnir, "Sample_ID", "Before", 1);
T_swir = movevars(T_swir, "Sample_ID", "Before", 1);

%% Merge by common sample ID
T_merged = innerjoin(T_vnir, T_swir, "Keys", "Sample_ID");

%% Check result
size(T_merged)
head(T_merged)

mergedFile = fullfile(outputDir, 'Merged_VNIR_SWIR_Summary.mat');
save(mergedFile, 'T_merged');

%% Creating merged cells

%% Create combined VNIR + SWIR spectra column
% Rule:
% Combined spectrum = full VNIR spectrum + SWIR spectrum from index 9 onwards

clear; clc;

%% Load merged table

fileIn = fullfile(outputDir, 'Merged_VNIR_SWIR_Summary.mat');
S = load(fileIn, 'T_merged');
T = S.T_merged;

%% Check required columns exist
requiredVars = ["VNIR_Mean_Spectrum", "SWIR_Mean_Spectrum"];

if ~all(ismember(requiredVars, T.Properties.VariableNames))
    error("The table must contain VNIR_Mean_Spectrum and SWIR_Mean_Spectrum columns.");
end

%% Create new columns
nSamples = height(T);

Combined_Spectrum = cell(nSamples, 1);
Number_Of_Combined_Wavelengths = nan(nSamples, 1);

for i = 1:nSamples

    % Extract spectra for sample i
    vnirSpectrum = T.VNIR_Mean_Spectrum{i};
    swirSpectrum = T.SWIR_Mean_Spectrum{i};

    % Force both spectra to be row vectors
    vnirSpectrum = vnirSpectrum(:).';
    swirSpectrum = swirSpectrum(:).';

    % Combine full VNIR with SWIR from 9th index onwards
    combinedSpectrum = [vnirSpectrum, swirSpectrum(9:end)];

    % Store combined spectrum and number of wavelengths
    Combined_Spectrum{i} = combinedSpectrum;
    Number_Of_Combined_Wavelengths(i) = numel(combinedSpectrum);

end

%% Add new columns to table
T.Combined_Spectrum = Combined_Spectrum;
T.Number_Of_Combined_Wavelengths = Number_Of_Combined_Wavelengths;

%% Check result
disp(T(:, ["Sample_ID", ...
           "VNIR_Number_Of_Wavelengths", ...
           "SWIR_Number_Of_Wavelengths", ...
           "Number_Of_Combined_Wavelengths"]));

%% Save updated table
T_merged = T;
save('Merged_VNIR_SWIR_Summary_With_Combined.mat', 'T_merged');
