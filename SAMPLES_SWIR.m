cd('C:\Users\ssomani\Desktop\Dinmukhamet\Master-Thesis-Automating-raw-data-sampling\Analysis')
clear variables 
close all 

%%  Loading the data

% Getting the source folder
srcFolder = fullfile(getenv('USERPROFILE'), 'Desktop','Dinmukhamet', 'Data');

% Finding all .mat files
allFiles = dir(fullfile(srcFolder, '*.mat'));

fileNames = {allFiles.name};


%Keeping the SWIR  files only 
isSWIR = contains(fileNames, '_SWIR', 'IgnoreCase', true);
% Excluding the MIX files
isMix = contains(fileNames, '_MIX_', 'IgnoreCase', true);

swirFiles = allFiles(isSWIR & ~isMix);

outputDir = fullfile(pwd, 'analysis_outputs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

SWIR_Summary = table( ...
    cell(0,1), ...
    cell(0,1), ...
    cell(0,1), ...
    zeros(0,1), ...
    'VariableNames', {'Sample_Name', 'Mean_Spectrum', 'Wavelengths', 'Number_Of_Wavelengths'} ...
);



%% Looping over 
for f = 1:numel(swirFiles)
    filePath = fullfile(swirFiles(f).folder, swirFiles(f).name);
    loadedData = load(filePath);
    
    
    im_reflectance = loadedData.im_reflectance;
    WL = loadedData.WL;


    mean_image = mean(im_reflectance, 3);
    figure, imshow(mean_image, [])

    [~, cubeNameRaw, ~] = fileparts(swirFiles(f).name);
    cubeName = regexprep(cubeNameRaw, '_SWIR_.*$', '_SWIR');
    cubeLabel = replace(cubeName, "_", " ");
    title("Mean Image Across All Wavelength Bands for SWIR: " + cubeLabel, 'Interpreter', 'none'); 

    im_uf = unfold(im_reflectance);
    im_uf_sd = std(im_uf, 0, 2);
    im_sd = reshape(im_uf_sd, size(im_reflectance, 1), size(im_reflectance, 2));

    figure
    subplot(1,2,1)
    imshow(im_sd,[]);

    subplot(1,2,2)
    histogram(im_sd,100);
    title('Pixel Intensities of standard deviation')
    ylim([0 8e4]);
    ylabel('Number of Pixels');
    xlabel('Standard Deviation of Pixel Intensities Across All Wavelength Bands');
    title('Pixel Intensities of standard deviation')


    mask_mix = im_sd > 2.3;

    figure, imshow(mask_mix, [])
    title('Mask of Pixels with Standard Deviation > 2.3')

    figure, imshow(mask_mix .* mean_image, [])
    title('Mean Image with Mask Applied')

    mask2 = imfill(mask_mix, "holes");
    figure, imshow(mask2,[])

    figure, imshow(mask2,[])
    title('Mask after Filling Holes')
    
    figure, imshow(mask2 .* mean_image, [])
    title('Mean Image with Mask Applied')

    L = bwlabel(mask2);
    figure, imshow(L,[])
    colormap('jet')

        x = zeros(1,20);
    y = zeros(1,20);

    for i = 1:20
        fprintf('Select the central region of sample %d\n', i);
        [xi, yi, ~] = impixel(im_reflectance(:,:,10),[]);
        x(i) = round(xi(1));
        y(i) = round(yi(1));
    end

    L_new = zeros(size(L));

    for i = 1:20
        L_old = L(y(i), x(i));
        L_new(L == L_old) = i;
    end

    

    % Finding mean spectra for each sample and creating a summary table
  nSamples = max(L_new(:));
  nBands   = size(im_reflectance, 3);

  im_reshaped = reshape(im_reflectance, [], nBands);
  meanSpectra = nan(nSamples, nBands);
    samplePixelCoords = cell(nSamples, 1);

  for i = 1:nSamples
        sampleMask = (L_new == i);
        samplePixelCoords{i} = find(sampleMask);

        samplePixels = im_reshaped(sampleMask(:), :);
    samplePixels = samplePixels(any(samplePixels, 2), :);
    if ~isempty(samplePixels)
        meanSpectra(i, :) = mean(samplePixels, 1);
    end
 end
  sampleNames = cubeName + "_sample_" + compose("%02d", (1:nSamples)');
 
  figure, imshow(L_new > 0, [])
    colormap('jet');
        title("New Label Matrix with 20 Samples Labeled for " + cubeLabel, 'Interpreter', 'none');

 meanSpectrumCell = mat2cell(meanSpectra, ones(nSamples,1), nBands);

 wavelengthCell = repmat({WL}, nSamples, 1);

 spectraSummaryTable = table( ...
    cellstr(sampleNames), meanSpectrumCell,wavelengthCell, ...
    repmat(nBands, nSamples, 1), 'VariableNames', {'Sample_Name', 'Mean_Spectrum', 'Wavelengths', 'Number_Of_Wavelengths'}); %#ok<NOPTS>

 SWIR_Summary = [SWIR_Summary; spectraSummaryTable];

end

save(fullfile(outputDir, 'SWIR_Summary.mat'), 'SWIR_Summary');





























