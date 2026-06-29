clear variables
close all

cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Analysis')


%% Load sample mean spectra from MIX_VNIR and MIX_SWIR

vnirData = load('MIX_VNIR_sampleMean.mat');
swirData = load('MIX_SWIR_sampleMean.mat');

sampleMean_VNIR = vnirData.sampleMean_VNIR;
sampleMean_SWIR = swirData.sampleMean_SWIR;

WL_vnir = vnirData.WL(:);
WL_swir = swirData.WL(:);
WL_combined = [WL_vnir; WL_swir(9:end)];

%% Building the table for classification

nSamples = numel(sampleMean_VNIR);

X_vnir_mean = zeros(nSamples, numel(WL_vnir));
X_swir_mean = zeros(nSamples, numel(WL_swir));
Y = strings(nSamples, 1);



%% Loop through samples to build data matrices and combined sample structure
MIX_sampleMean = struct();

for sampleID = 1:nSamples
    X_vnir_mean(sampleID,:) = sampleMean_VNIR(sampleID).meanSpectrum;
    X_swir_mean(sampleID,:) = sampleMean_SWIR(sampleID).meanSpectrum;
    Y(sampleID) = sampleMean_VNIR(sampleID).className;

    MIX_sampleMean(sampleID).sampleID = sampleID;
    MIX_sampleMean(sampleID).sampleName = sampleMean_VNIR(sampleID).sampleName;
    MIX_sampleMean(sampleID).classID = sampleMean_VNIR(sampleID).classID;
    MIX_sampleMean(sampleID).className = sampleMean_VNIR(sampleID).className;
    MIX_sampleMean(sampleID).wavelengths_VNIR = WL_vnir;
    MIX_sampleMean(sampleID).absorbance_VNIR = X_vnir_mean(sampleID,:);
    MIX_sampleMean(sampleID).wavelengths_SWIR = WL_swir;
    MIX_sampleMean(sampleID).absorbance_SWIR = X_swir_mean(sampleID,:);
    MIX_sampleMean(sampleID).wavelengths = WL_combined;
    MIX_sampleMean(sampleID).meanSpectrum = [X_vnir_mean(sampleID,:) X_swir_mean(sampleID,9:end)];
    MIX_sampleMean(sampleID).unfoldedSpectra_VNIR = sampleMean_VNIR(sampleID).unfoldedSpectra;
    MIX_sampleMean(sampleID).unfoldedSpectra_SWIR = sampleMean_SWIR(sampleID).unfoldedSpectra;
    MIX_sampleMean(sampleID).nPixels_VNIR = sampleMean_VNIR(sampleID).nPixels;
    MIX_sampleMean(sampleID).nPixels_SWIR = sampleMean_SWIR(sampleID).nPixels;
end





%% Save prepared data

save('MIX_classification_inputs.mat', ...
    'MIX_sampleMean')

disp("Saved MIX_classification_inputs.mat");



