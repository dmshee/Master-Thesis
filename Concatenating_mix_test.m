% cd('')
clear variables 
close all 

% scriptDir = "";

%% Load the VNIR and SWIR summary tables
load(fullfile(scriptDir,'MIX_VNIR_Summary.mat'));
load(fullfile(scriptDir,'MIX_SWIR_Summary.mat'));

A = MIX_VNIR_Summary;
B = MIX_SWIR_Summary;

A.Properties.VariableNames{'Mean_Spectrum'} = 'VNIR_Mean_Spectrum';
A.Properties.VariableNames{'Wavelengths'} = 'VNIR_Bands';

B.Properties.VariableNames{'Mean_Spectrum'} = 'SWIR_Mean_Spectrum';
B.Properties.VariableNames{'Wavelengths'} = 'SWIR_Bands';


A.Sample_Key = regexprep(string(A.Sample_Name), '_VNIR_', '_');
B.Sample_Key = regexprep(string(B.Sample_Name), '_SWIR_', '_');


assert(height(A) == numel(unique(A.Sample_Key)), 'A has duplicate sample keys');
assert(height(B) == numel(unique(B.Sample_Key)), 'B has duplicate sample keys');

% One-to-one join
M = innerjoin(A, B, 'Keys', 'Sample_Key');

overlapBands = 8;
assert(all(cellfun(@numel, M.SWIR_Mean_Spectrum) > overlapBands), ...
	'Some SWIR spectra are shorter than overlapBands.');
assert(all(cellfun(@numel, M.VNIR_Mean_Spectrum) > overlapBands), ...
	'Some VNIR spectra are shorter than overlapBands.');

M.Combined_Mean_Spectrum = cell(height(M), 1);
M.Scale_Factor = nan(height(M), 1);
M.Mean_Shift = nan(height(M), 1);

for i = 1:height(M)
	v = M.VNIR_Mean_Spectrum{i};
	s = M.SWIR_Mean_Spectrum{i};

	vOverlap = v(end-overlapBands+1:end);
	sOverlap = s(1:overlapBands);

	vMean = mean(vOverlap);
	sMean = mean(sOverlap);
	vStd = std(vOverlap);
	sStd = std(sOverlap);

	if sStd <= eps
		scaleFactor = 1;
	else
		scaleFactor = vStd / sStd;
	end

	sOverlapAutoscaled = (sOverlap - sMean) * scaleFactor + vMean;
	overlapBlend = (vOverlap + sOverlapAutoscaled) / 2;
	M.Combined_Mean_Spectrum{i} = [v(1:end-overlapBands), overlapBlend, s(overlapBands+1:end)];

	M.Scale_Factor(i) = scaleFactor;
	M.Mean_Shift(i) = vMean - sMean * scaleFactor;
end

M.Number_Of_Wavelengths = cellfun(@numel, M.Combined_Mean_Spectrum);

Tbl_combined_mix = M(:, {'Sample_Name_A', 'Sample_Name_B', 'Combined_Mean_Spectrum'});

%% Adding the labels
sampleName = string(Tbl_combined_mix.Sample_Name_A);
sampleType = regexprep(sampleName, '_sample_\d+$', '');
nameParts = split(sampleName, '_');
sampleClass = nameParts(:,2);
[sampleTypes, ~, grp] = unique(sampleClass, 'stable');
Tbl_combined_mix.sampleType = sampleType;
Tbl_combined_mix.sampleClass = sampleClass;

save(fullfile(scriptDir, 'Tbl_combined_mix.mat'), 'Tbl_combined_mix', '-v7.3');



%% Plotting the full combined spectra to ensure if they are consistent 
allLengths = cellfun(@numel, Tbl_combined_mix.Combined_Mean_Spectrum);
assert(isscalar(unique(allLengths)), 'Combined spectra lengths are inconsistent across samples.');

combinedSpectra = vertcat(Tbl_combined_mix.Combined_Mean_Spectrum{:});
bandIndex = 1:size(combinedSpectra, 2);

figure('Color', 'w', 'Name', 'Full Combined Spectra (Mix)');
plot(bandIndex, combinedSpectra', 'LineWidth', 1);
hold on;
hMean = plot(bandIndex, mean(combinedSpectra, 1), 'k--', 'LineWidth', 4);
hold off;
xlabel('Band Index');
ylabel('Reflectance / Intensity');
title('Full Combined VNIR-SWIR Spectra for Mix Test Samples');
grid on;

figuresDir = 'C:\Users\ssomani\Desktop\Dinmukhamet\Master-Thesis-Automating-raw-data-sampling';
outFig = fullfile(figuresDir, 'Figures', 'Full_Combined_Spectra_Mix_Test.png');
exportgraphics(gcf, outFig, 'Resolution', 300);
