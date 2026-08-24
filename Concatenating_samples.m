% cd('')
clear variables 
close all 


% sourcedir = "";

srcFolder = uigetdir(sourcedir,'source folder');
files = dir(fullfile(srcFolder, '*.mat'));

for i = 1:numel(files)
    filePath = fullfile(srcFolder, files(i).name);
    load(filePath);
end 


%% Concatenating the VNIR and SWIR data
A = VNIR_Summary;
B = SWIR_Summary;

A.Properties.VariableNames{'Mean_Spectrum'} = 'VNIR_Mean_Spectrum';
A.Properties.VariableNames{'Wavelengths'} = 'VNIR_Bands';


B.Properties.VariableNames{'Mean_Spectrum'} = 'SWIR_Mean_Spectrum';
B.Properties.VariableNames{'Wavelengths'} = 'SWIR_Bands';



A.Sample_Key = regexprep(string(A.Sample_Name), '_VNIR_', '_');
B.Sample_Key = regexprep(string(B.Sample_Name), '_SWIR_', '_');

assert(height(A) == numel(unique(A.Sample_Key)), 'A has duplicate sample keys');
assert(height(B) == numel(unique(B.Sample_Key)), 'B has duplicate sample keys');


M = innerjoin(A, B, 'Keys', 'Sample_Key');

overlapBands = 7;

M.Combined_Mean_Spectrum = cell(height(M), 1);
M.Combined_Wavelengths = cell(height(M), 1);

for i = 1:height(M)
    vnirSpectrum = M.VNIR_Mean_Spectrum{i};
    swirSpectrum = M.SWIR_Mean_Spectrum{i};

    vnirOverlap = vnirSpectrum(end-overlapBands+1:end);
    swirOverlap = swirSpectrum(1:overlapBands);

    vnirMean = mean(vnirOverlap);
    swirMean = mean(swirOverlap);
    vnirStd = std(vnirOverlap);
    swirStd = std(swirOverlap);

    if swirStd <= eps
        scaleFactor = 1;
    else
        scaleFactor = vnirStd / swirStd;
    end

    swirOverlapAutoscaled = (swirOverlap - swirMean) * scaleFactor + vnirMean;
    overlapBlend = (vnirOverlap + swirOverlapAutoscaled) / 2;

    M.Combined_Mean_Spectrum{i} = [vnirSpectrum(1:end-overlapBands), overlapBlend, swirSpectrum(overlapBands+1:end)];

    M.Scale_Factor(i) = scaleFactor;
    M.Mean_Shift(i) = vnirMean - swirMean * scaleFactor;

    %% Creating combined wavelengths indexes
    vW = M.VNIR_Bands{i};
    sW = M.SWIR_Bands{i};

    vWOverlap = vW(end-overlapBands+1:end);
    sWOverlap = sW(1:overlapBands);

    overlapWBlend = (vWOverlap + sWOverlap) / 2;
    M.Combined_Wavelengths{i} = [vW(1:end-overlapBands), overlapWBlend, sW(overlapBands+1:end)];
end

M.Combined_Number_Of_Wavelengths = cellfun(@numel, M.Combined_Mean_Spectrum);

Tbl_combined_samples = M(:, {'Sample_Name_A', 'Sample_Name_B', 'Combined_Mean_Spectrum', 'Combined_Wavelengths'});
save(fullfile(sourcedir, 'Tbl_combined_samples.mat'), 'Tbl_combined_samples', '-v7.3');

%% Plotting the full combined spectra
combinedSpectra = vertcat(Tbl_combined_samples.Combined_Mean_Spectrum{:});
bandIndex = Tbl_combined_samples.Combined_Wavelengths{1};


figure('Color', 'w', 'Name', 'Full Combined Spectra');
plot(bandIndex, combinedSpectra', 'LineWidth', 1);
hold on;
hMean = plot(bandIndex, mean(combinedSpectra, 1),'k--', 'LineWidth', 4);
hold off;
xlabel('Wavelength (nm)');
ylabel('Reflectance (a.u.)');
title('Full Combined VNIR-SWIR Spectra for Samples (vertical concatenation)');
grid on;

exportgraphics(gcf, fullfile(sourcedir, 'Full_Combined_Spectra.png'), 'Resolution', 300);


