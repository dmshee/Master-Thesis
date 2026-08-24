% cd('')
clear variables
close all


%% Saving the information 
% srcFolder = fullfile(getenv('USERPROFILE'),'Desktop','','Data');
allFiles  = dir(fullfile(srcFolder, '*.mat'));
fileNames = {allFiles.name};

isSWIR = contains(fileNames, '_SWIR', 'IgnoreCase', true);
isMix  = contains(fileNames, '_MIX_',  'IgnoreCase', true);
swirFiles = allFiles(isSWIR & isMix);

outputDir = pwd;


summaryPath = fullfile(outputDir, 'MIX_SWIR_Summary.mat');


MIX_SWIR_Summary = table( ...
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

  classNames   = {'HLTHY', 'BEM','BEN','PEN'};
  classSamples = [4,7,10,13,17; 3,6,9,12,15; 8,11,14,18,20; 1,2,5,16,19];

  sample_to_class = zeros(1,20);
   for classID = 1:numel(classNames)
     sample_to_class(classSamples(classID,:)) = classID;
   end
  
     %% Printing the mean image
    mean_image = mean(im_reflectance, 3);
    figure, imshow(mean_image, [])

    [~, cubeNameRaw, ~] = fileparts(swirFiles(f).name);
    cubeName = string(regexprep(cubeNameRaw, '_SWIR_.*$', '_SWIR'));
    cubeLabel = replace(cubeName, "_", " ");
    title("Mean Image Across All Wavelength Bands for SWIR: " + cubeLabel, 'Interpreter', 'none'); 

    %% Unfolding the spectra
    Xuf   = unfold(im_reflectance);
    im_sd = reshape(std(Xuf,0,2), size(im_reflectance,1), size(im_reflectance,2));

    % displaying the histogram of standard deviation 
    figure, histogram(im_sd, 100);
    ylim([0 5e4]);
    ylabel('Number of Pixels');
    xlabel('Standard Deviation of Pixel Intensities Across All Wavelength Bands');
    title('Pixel Intensities of standard deviation')

    

    
    mask1 = im_sd > 2.25;
    mask2 = imfill(mask1, "holes");

    figure;
    imshow(mask2, []);
    title("Mask after applying imfill for SWIR: " + cubeLabel, 'Interpreter', 'none');
    drawnow;

    % displaying the masks together in one figure
    figure;
    imshow(mask1, []);
    title("Mask of Pixels with Standard Deviation > 2.25 for SWIR: " + cubeLabel, 'Interpreter', 'none');



    L     = bwlabel(mask2);
    
    %% Selecting the samples
    x = zeros(1,20);
    y = zeros(1,20);
    for i = 1:20
        fprintf('File %d/%d (%s): Select center of sample %d.\n', f, numel(swirFiles), swirFiles(f).name, i);
        [xi, yi, ~] = impixel(im_reflectance(:,:,10), []);
        x(i) = round(xi(1));
        y(i) = round(yi(1));
    end

    L_new = zeros(size(L));
    for i = 1:20
        L_old = L(y(i), x(i));
        L_new(L == L_old) = i;
    end
    
 

   %% Displaying the class map
    class_map = zeros(size(L_new));
   for s = 1:20
    class_map(L_new == s) = sample_to_class(s);
   end

   figure,
   imshow(class_map, []);
   colormap('jet');
   colorbar('Ticks', 1:length(classNames), 'TickLabels', classNames);
   title("Class Map of Samples for SWIR: " + cubeLabel, 'Interpreter', 'none'); 
 

   %% Finding the mean spectra for each sample and saving to summary table
   
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
  
  rawName = string(cubeNameRaw);   % original file stem (before .mat)

  % Detect position (up/down) from cube name
  position = "unknown";
  if contains(lower(rawName), "_up")
      position = "up";
  elseif contains(lower(rawName), "_down")
      position = "down";
  end

  % Map name = part before "_up" or "_down" 
  mapName = regexprep(rawName, '_(up|down).*$', '', 'ignorecase');
  if strlength(mapName) == 0
      mapName = cubeName;
  end

  
  sampleIDs = (1:nSamples)';
  classIDs = sample_to_class(sampleIDs);
  sampleClassNames = string(classNames(classIDs))';
  mapPrefixWithClass = strings(nSamples,1);
  for i = 1:nSamples
      mapPrefixWithClass(i) = regexprep(mapName, '(?i)MIX', sampleClassNames(i));
  end

  
  sampleNames = mapPrefixWithClass + "_" + position + "_sample_" + compose("%02d", sampleIDs);

  figure, imshow(L_new > 0, [])
    colormap('jet');
        title("New Label Matrix with 20 Samples Labeled for " + cubeLabel, 'Interpreter', 'none');


 meanSpectrumCell = mat2cell(meanSpectra, ones(nSamples,1), nBands);

 wavelengthCell = repmat({WL}, nSamples, 1);

  
 spectraSummaryTable = table( ...
    cellstr(sampleNames), meanSpectrumCell, wavelengthCell, ...
    repmat(nBands, nSamples, 1), ...
    'VariableNames', {'Sample_Name', 'Mean_Spectrum', 'Wavelengths', 'Number_Of_Wavelengths'} ...
 );

 MIX_SWIR_Summary = [MIX_SWIR_Summary; spectraSummaryTable]; %#ok<AGROW>

end

save(summaryPath, 'MIX_SWIR_Summary', '-v7.3');
fprintf('Saved %d rows to %s\n', height(MIX_SWIR_Summary), summaryPath);




%% Showing the map 
%Building class map image from sample labels
class_map = zeros(size(L_new));
for s = 1:20
    class_map(L_new == s) = sample_to_class(s);
end

% Displaying the class map
figure,
imshow(class_map, []);
colormap('jet');
colorbar('Ticks', 1:length(classNames), 'TickLabels', classNames);
title("Class Map of Samples for SWIR: " + cubeLabel, 'Interpreter', 'none'); 