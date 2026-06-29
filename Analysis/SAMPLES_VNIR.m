cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Analysis') 
close all 

load('BEN_VNIR.mat')

mean_image = mean(im_reflectance, 3);
figure, imshow(mean_image, [])
title('Mean Image Across All Wavelength Bands')

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
title('Pixel Intensities of sd')

% Applying the mask to the mean image
mask_mix = im_sd>2.916;
figure, imshow(mask_mix.*mean_image,[])
title('Mean Image with Mask Applied')


figure, imshow(mask_mix.*mean_image,[])
title('Mean Image with Final Mask Applied')

%Filling in the holes in the mask 
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
    [ xi, yi, ~] = impixel(im_reflectance(:,:,150),[]);
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

sampleNames = "BEN_VNIR_Sample_" + compose("%02d", (1:nSamples)');

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

fileOut = fullfile(outputDir, 'mean_sample_summary.mat');

% New batch computed
newTable = table( ...
    cellstr(sampleNames), ...
    meanSpectrumCell, ...
    repmat(nBands, nSamples, 1), ...
    'VariableNames', {'Sample_Name', 'Mean_Spectrum', 'Number_Of_Wavelengths'} ...
);

% Append if file exists, otherwise create
if isfile(fileOut)
    S = load(fileOut, 'mean_sample_summary');
    oldTable = S.mean_sample_summary;

    % Optional: avoid duplicate sample names
    keep = ~ismember(newTable.Sample_Name, oldTable.Sample_Name);
    newTable = newTable(keep, :);

    mean_sample_summary = [oldTable; newTable];
else
    mean_sample_summary = newTable;
end

save(fileOut, 'mean_sample_summary');

