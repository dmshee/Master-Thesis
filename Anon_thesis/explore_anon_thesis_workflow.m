%% =========================================================================
%% ANONYMIZED THESIS WORKFLOW: EXPLORATORY ANALYSIS SCRIPT
%% =========================================================================
%% 
%% This script performs hyperspectral image analysis on kidney bean samples:
%% - Preprocessing and segmentation
%% - PCA comparison across preprocessing methods
%% - Unsupervised clustering (k-means)
%% - Supervised classification (PLS-DA) for age and germination
%%
%% Run this script section by section to inspect intermediate data and figures.
%% All major intermediate variables remain visible in the workspace.

%% =========================================================================
%% 1. Setup and file paths
%% =========================================================================

clearvars;
close all;
clc;

fprintf('\n========== ANONYMIZED THESIS WORKFLOW: EXPLORATORY MODE ==========\n\n');

baseDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(baseDir);
addpath(baseDir);
addpath(parentDir);

% Create output directory if it does not exist
outputDir = fullfile(baseDir, 'analysis_outputs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% =========================================================================
%% 2. Exploratory Parameters
%% =========================================================================
%% These parameters control the analysis workflow. Modify them to explore
%% different segmentation, preprocessing, or classification strategies.

% RGB visualization bands (nanometers)
rgbWavelengths = [650, 550, 450];

% Ratio mask bands (nanometers)
ratioHighWavelength = 900;
ratioLowWavelength = 650;

% First-pass PCA: number of components to compute and visualize
nFirstPassPCs = 3;

% PCA-guided mask refinement: which components to use
pcForSecondMask = 2;
pcForThirdMask = 3;

% Mask refinement support region: pixels in [centralRegionFraction(1)*H : centralRegionFraction(2)*H]
centralRegionFraction = [0.3, 0.7];

% Score-based mask thresholds: minimum number of valid pixels
minimumScorePixels = 10;
minimumRetainedPixelsAbsolute = 50;
minimumRetainedPixelsFraction = 0.2;

% K-means clustering
nKmeansClusters = 2;
kmeansMaxIterations = 50;
randomSeed = 7;

% Classification data split
calibrationPerClass = 75;
testPerClass = 25;

% PLS-DA parameters
maxPlsComponentsCap = 20;
plsdaThreshold = 0.5;

% Figure output control
saveFigures = true;
closeFiguresAfterSaving = false;
saveResults = true;

fprintf('Parameters loaded.\n\n');

%% =========================================================================
%% 3. Load VNIR image cube
%% =========================================================================
cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Anon_thesis')
cubeFilePath = fullfile(baseDir, 'nb_s6_vnir_2.mat');
fprintf('Loading VNIR cube from: %s\n', cubeFilePath);

cubeData = load(cubeFilePath);
if ~isfield(cubeData, 'im')
    error('Expected variable "im" in nb_s6_vnir_2.mat, but it was not found.');
end

cube = double(cubeData.im);
fprintf('  Cube dimensions: %d x %d x %d (height x width x bands)\n', ...
    size(cube, 1), size(cube, 2), size(cube, 3));

%% =========================================================================
%% 4. Load bean-level ground-truth data
%% =========================================================================

gtFilePath = fullfile(baseDir, 'kidney_bean_data_GT.mat');
fprintf('Loading ground-truth data from: %s\n', gtFilePath);

gtData = load(gtFilePath);
if ~isfield(gtData, 'kidney_bean_data')
    error('Expected variable "kidney_bean_data" in kidney_bean_data_GT.mat, but it was not found.');
end

beanData = gtData.kidney_bean_data;
beanData = beanData(:)';  % Ensure row vector
fprintf('  Number of beans: %d\n', numel(beanData));

%% =========================================================================
%% 5. Build bean summary table
%% =========================================================================
%% Extract metadata from bean data structure: name, age (new/old), germination.

summaryTable = build_summary_table(beanData);

fprintf('\nBean summary:\n');
fprintf('  New beans: %d\n', sum(summaryTable.IsNew));
fprintf('  Old beans: %d\n', sum(~summaryTable.IsNew));
fprintf('  Germinated: %d\n', sum(summaryTable.Germination == 1));
fprintf('  Non-germinated: %d\n', sum(summaryTable.Germination == 0));

disp(head(summaryTable, 5));

%% =========================================================================
%% 6. Extract wavelength vectors
%% =========================================================================

vnirWavelengths = beanData(1).VNIR_WL(:);
fullWavelengths = beanData(1).Full_WL(:);

fprintf('\nWavelength vectors:\n');
fprintf('  VNIR bands: %d (range: %.0f–%.0f nm)\n', ...
    numel(vnirWavelengths), vnirWavelengths(1), vnirWavelengths(end));
fprintf('  Full-range bands: %d (range: %.0f–%.0f nm)\n', ...
    numel(fullWavelengths), fullWavelengths(1), fullWavelengths(end));

%% =========================================================================
%% 7. Initial RGB / pseudo-RGB visualization
%% =========================================================================
%% Create a pseudo-RGB image using three bands closest to target wavelengths.
%% This provides a familiar visual reference for the hyperspectral data.

redIdx = nearest_band(vnirWavelengths, rgbWavelengths(1));
greenIdx = nearest_band(vnirWavelengths, rgbWavelengths(2));
blueIdx = nearest_band(vnirWavelengths, rgbWavelengths(3));

rgbImage = cat(3, cube(:, :, redIdx), cube(:, :, greenIdx), cube(:, :, blueIdx));
rgbImage = scale_01(rgbImage);

figure('Name', 'Initial RGB visualization', 'Color', 'w');
imshow(rgbImage);
title('Pseudo-RGB composite (650 nm, 550 nm, 450 nm)');

maybe_save_figure(gcf, 'initial_rgb.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 8. First-pass 900/650 nm ratio mask
%% =========================================================================
%% The 900/650 nm ratio is used as a first-pass contrast image because it
%% separates bean tissue from background more clearly than a single band.
%% An Otsu threshold is applied to convert the continuous ratio image into
%% a binary mask.

fprintf('\nComputing first-pass ratio mask...\n');

highIdx = nearest_band(vnirWavelengths, ratioHighWavelength);
lowIdx = nearest_band(vnirWavelengths, ratioLowWavelength);

fprintf('  High band: %d nm (band index %d)\n', ratioHighWavelength, highIdx);
fprintf('  Low band: %d nm (band index %d)\n', ratioLowWavelength, lowIdx);

% Compute ratio image with numerical stability
ratioImage = cube(:, :, highIdx) ./ max(cube(:, :, lowIdx), eps);
ratioImage = scale_01(ratioImage);

% Apply Otsu threshold
ratioThreshold = otsu_threshold(ratioImage(:));
fprintf('  Otsu threshold: %.4f\n', ratioThreshold);

maskA = ratioImage >= ratioThreshold;
maskB = ratioImage < ratioThreshold;

% Choose polarity based on central support (beans are typically central)
centerRows = round(size(cube, 1) * centralRegionFraction(1)) : ...
             round(size(cube, 1) * centralRegionFraction(2));
centerCols = round(size(cube, 2) * centralRegionFraction(1)) : ...
             round(size(cube, 2) * centralRegionFraction(2));
centerSupport = false(size(maskA));
centerSupport(centerRows, centerCols) = true;

if nnz(maskA & centerSupport) >= nnz(maskB & centerSupport)
    firstMask = maskA;
    fprintf('  Chose maskA (ratio >= threshold)\n');
else
    firstMask = maskB;
    fprintf('  Chose maskB (ratio < threshold)\n');
end

fprintf('  First-pass mask: %d pixels\n', nnz(firstMask));

% Visualize ratio mask and first-pass mask
figure('Name', 'First-pass ratio mask', 'Color', 'w');

subplot(2, 2, 1);
imshow(rgbImage);
title('Pseudo-RGB composite');

subplot(2, 2, 2);
imagesc(ratioImage);
axis image off;
colorbar;
title(sprintf('Ratio image (%.0f / %.0f nm)', ratioHighWavelength, ratioLowWavelength));

subplot(2, 2, 3);
hold on;
histogram(ratioImage(:), 100, 'FaceAlpha', 0.7, 'EdgeAlpha', 0);
xline(ratioThreshold, 'r-', 'LineWidth', 2);
hold off;
xlabel('Ratio value');
ylabel('Count');
legend('Histogram', sprintf('Otsu threshold (%.4f)', ratioThreshold));
title('Ratio image histogram with Otsu threshold');

subplot(2, 2, 4);
imagesc(firstMask);
axis image off;
title(sprintf('First-pass mask (%d pixels)', nnz(firstMask)));

maybe_save_figure(gcf, 'mask_01_ratio_mask.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 9. First-pass PCA on masked spectra
%% =========================================================================
%% Extract spectra from the first-pass mask and perform PCA. This serves
%% two purposes: (1) dimensionality reduction for visualization, and
%% (2) identification of major spectral components to guide mask refinement.

fprintf('\nPerforming first-pass PCA on masked spectra...\n');

maskedSpectraFirst = extract_masked_spectra(cube, firstMask);
fprintf('  First-pass masked spectra: %d x %d (pixels x bands)\n', ...
    size(maskedSpectraFirst, 1), size(maskedSpectraFirst, 2));

[loadingsFirst, scoresFirst, eigenvaluesFirst] = PCA(maskedSpectraFirst);

explainedFirst = 100 * eigenvaluesFirst(:) / sum(eigenvaluesFirst);
fprintf('  Explained variance (first 5 PCs): %.2f, %.2f, %.2f, %.2f, %.2f %%\n', ...
    explainedFirst(1:min(5, numel(explainedFirst))));

% Reconstruct PCA score images
scoreImagesFirst = reconstruct_score_images(firstMask, scoresFirst, nFirstPassPCs);

% Visualize first-pass PCA
figure('Name', 'First-pass PCA', 'Color', 'w');

subplot(2, 4, 1);
plot(explainedFirst(1:min(15, numel(explainedFirst))), 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Principal component');
ylabel('Explained variance (%)');
title('Explained variance');

for pc = 1:3
    subplot(2, 4, pc + 1);
    plot(vnirWavelengths, loadingsFirst(:, pc), 'LineWidth', 1.2);
    grid on;
    xlabel('Wavelength (nm)');
    ylabel('Loading');
    title(sprintf('PC%d loading (%.1f%%)', pc, explainedFirst(pc)));
end

for pc = 1:3
    subplot(2, 4, pc + 4);
    imagesc(scoreImagesFirst(:, :, pc));
    axis image off;
    colorbar;
    title(sprintf('PC%d score image', pc));
end

maybe_save_figure(gcf, 'pca_01_first_pass.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 10. PCA-guided mask refinement using PC2 and PC3
%% =========================================================================
%% The PC2 and PC3 score images often reveal internal structure within the
%% bean that indicates region of interest. Secondary masks are derived using
%% Otsu thresholding on these score images, then combined with the first
%% mask to refine the region of interest progressively.

fprintf('\nRefining mask using PC%d and PC%d...\n', pcForSecondMask, pcForThirdMask);

% Second mask: use PC2 score image
pc2ScoreImage = scoreImagesFirst(:, :, pcForSecondMask);
[secondMask, pc2Threshold] = score_mask(...
    pc2ScoreImage, firstMask, ...
    centralRegionFraction, minimumScorePixels, ...
    minimumRetainedPixelsAbsolute, minimumRetainedPixelsFraction);

fprintf('  PC%d threshold: %.4f, second mask: %d pixels\n', ...
    pcForSecondMask, pc2Threshold, nnz(secondMask));

% Third mask: use PC3 score image, with second mask as support
pc3ScoreImage = scoreImagesFirst(:, :, pcForThirdMask);
support2 = firstMask & secondMask;
[thirdMask, pc3Threshold] = score_mask(...
    pc3ScoreImage, support2, ...
    centralRegionFraction, minimumScorePixels, ...
    minimumRetainedPixelsAbsolute, minimumRetainedPixelsFraction);

fprintf('  PC%d threshold: %.4f, third mask: %d pixels\n', ...
    pcForThirdMask, pc3Threshold, nnz(thirdMask));

% Final mask: intersection of all three
finalMask = firstMask & secondMask & thirdMask;
fprintf('  Final mask: %d pixels\n', nnz(finalMask));

% Visualize mask refinement
figure('Name', 'Mask refinement using PC scores', 'Color', 'w');

subplot(2, 3, 1);
imagesc(firstMask);
axis image off;
title(sprintf('First mask\n(%d pixels)', nnz(firstMask)));

subplot(2, 3, 2);
imagesc(pc2ScoreImage);
axis image off;
colorbar;
title(sprintf('PC%d score image', pcForSecondMask));

subplot(2, 3, 3);
imagesc(secondMask);
axis image off;
title(sprintf('Second mask\n(%d pixels)', nnz(secondMask)));

subplot(2, 3, 4);
imagesc(pc3ScoreImage);
axis image off;
colorbar;
title(sprintf('PC%d score image', pcForThirdMask));

subplot(2, 3, 5);
imagesc(thirdMask);
axis image off;
title(sprintf('Third mask\n(%d pixels)', nnz(thirdMask)));

subplot(2, 3, 6);
imagesc(finalMask);
axis image off;
title(sprintf('Final mask\n(%d pixels)', nnz(finalMask)));

maybe_save_figure(gcf, 'mask_02_refinement.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 11. Final masked spectra extraction
%% =========================================================================
%% Extract spectra from the final refined mask. These spectra will be
%% preprocessed and analyzed separately.

fprintf('\nExtracting final masked spectra...\n');

maskedSpectra = extract_masked_spectra(cube, finalMask);
fprintf('  Final masked spectra: %d x %d (pixels x bands)\n', ...
    size(maskedSpectra, 1), size(maskedSpectra, 2));

%% =========================================================================
%% 12. Spectral preprocessing: raw, SNV, MSC
%% =========================================================================
%% Three preprocessing variants are created:
%% - Raw: no preprocessing
%% - SNV: Standard Normal Variate (removes baseline offset and scaling)
%% - MSC: Multiplicative Signal Correction (corrects for multiplicative effects)

fprintf('\nApplying spectral preprocessing...\n');

snvSpectra = SNV(maskedSpectra);
fprintf('  SNV spectra computed.\n');

mscSpectra = msc(maskedSpectra);
fprintf('  MSC spectra computed.\n');

% Compare preprocessing with example spectra
nExamplesPerGroup = 50;

exampleIndicesRaw = randperm(size(maskedSpectra, 1), min(nExamplesPerGroup, size(maskedSpectra, 1)));
exampleRaw = maskedSpectra(exampleIndicesRaw, :);
exampleSnv = snvSpectra(exampleIndicesRaw, :);
exampleMsc = mscSpectra(exampleIndicesRaw, :);

figure('Name', 'Preprocessing comparison', 'Color', 'w');

subplot(1, 3, 1);
plot(vnirWavelengths, exampleRaw', 'Color', [0.7, 0.7, 0.9], 'LineWidth', 0.5);
hold on;
plot(vnirWavelengths, mean(exampleRaw, 1), 'b', 'LineWidth', 2);
hold off;
grid on;
xlabel('Wavelength (nm)');
ylabel('Reflectance');
title('Raw spectra');

subplot(1, 3, 2);
plot(vnirWavelengths, exampleSnv', 'Color', [0.9, 0.7, 0.7], 'LineWidth', 0.5);
hold on;
plot(vnirWavelengths, mean(exampleSnv, 1), 'r', 'LineWidth', 2);
hold off;
grid on;
xlabel('Wavelength (nm)');
ylabel('SNV');
title('SNV-preprocessed spectra');

subplot(1, 3, 3);
plot(vnirWavelengths, exampleMsc', 'Color', [0.7, 0.9, 0.7], 'LineWidth', 0.5);
hold on;
plot(vnirWavelengths, mean(exampleMsc, 1), 'g', 'LineWidth', 2);
hold off;
grid on;
xlabel('Wavelength (nm)');
ylabel('MSC');
title('MSC-preprocessed spectra');

maybe_save_figure(gcf, 'preprocessing_comparison.png', saveFigures, closeFiguresAfterSaving, outputDir);

fprintf('  Preprocessing comparison figure shown.\n');

%% =========================================================================
%% 13. PCA comparison across preprocessing methods
%% =========================================================================
%% PCA is performed independently on raw, SNV, and MSC data. This allows
%% comparison of how preprocessing affects the principal component structure.

fprintf('\nPerforming PCA on all three preprocessing variants...\n');

[loadingsRaw, scoresRaw, eigenvaluesRaw] = PCA(maskedSpectra);
fprintf('  PCA on raw spectra: %d components\n', numel(eigenvaluesRaw));

[loadingsSnv, scoresSnv, eigenvaluesSnv] = PCA(snvSpectra);
fprintf('  PCA on SNV spectra: %d components\n', numel(eigenvaluesSnv));

[loadingsMsc, scoresMsc, eigenvaluesMsc] = PCA(mscSpectra);
fprintf('  PCA on MSC spectra: %d components\n', numel(eigenvaluesMsc));

% Reconstruct score images
scoreImagesRaw = reconstruct_score_images(finalMask, scoresRaw, 3);
scoreImagesSnv = reconstruct_score_images(finalMask, scoresSnv, 3);
scoreImagesMsc = reconstruct_score_images(finalMask, scoresMsc, 3);

% Plot PCA comparison for each preprocessing method
for prep_idx = 1:3
    if prep_idx == 1
        eigenvalues = eigenvaluesRaw;
        loadings = loadingsRaw;
        scoreImages = scoreImagesRaw;
        tag = 'raw';
    elseif prep_idx == 2
        eigenvalues = eigenvaluesSnv;
        loadings = loadingsSnv;
        scoreImages = scoreImagesSnv;
        tag = 'snv';
    else
        eigenvalues = eigenvaluesMsc;
        loadings = loadingsMsc;
        scoreImages = scoreImagesMsc;
        tag = 'msc';
    end
    
    explained = 100 * eigenvalues(:) / sum(eigenvalues);
    
    figure('Name', sprintf('PCA (%s)', tag), 'Color', 'w');
    
    subplot(2, 4, 1);
    plot(explained(1:min(15, numel(explained))), 'o-', 'LineWidth', 1.5);
    grid on;
    xlabel('Principal component');
    ylabel('Explained variance (%)');
    title('Explained variance');
    
    for pc = 1:3
        subplot(2, 4, pc + 1);
        plot(vnirWavelengths, loadings(:, pc), 'LineWidth', 1.2);
        grid on;
        xlabel('Wavelength (nm)');
        ylabel('Loading');
        title(sprintf('PC%d (%.1f%%)', pc, explained(pc)));
    end
    
    for pc = 1:2
        subplot(2, 4, pc + 4);
        imagesc(scoreImages(:, :, pc));
        axis image off;
        colorbar;
        title(sprintf('PC%d score', pc));
    end
    
    maybe_save_figure(gcf, sprintf('pca_%s.png', tag), saveFigures, closeFiguresAfterSaving, outputDir);
end

fprintf('  PCA comparison figures shown.\n');

%% =========================================================================
%% 14. K-means clustering comparison
%% =========================================================================
%% Unsupervised k-means clustering is applied to each preprocessing variant.
%% This reveals whether preprocessing affects the natural grouping within
%% the bean tissue.

fprintf('\nPerforming k-means clustering (%d clusters)...\n', nKmeansClusters);

rng(randomSeed);

clusterLabelsRaw = kmeans_simple(maskedSpectra, nKmeansClusters, kmeansMaxIterations);
clusterRaw = reconstruct_label_image(finalMask, clusterLabelsRaw);
fprintf('  Raw clustering complete.\n');

clusterLabelsSnv = kmeans_simple(snvSpectra, nKmeansClusters, kmeansMaxIterations);
clusterSnv = reconstruct_label_image(finalMask, clusterLabelsSnv);
fprintf('  SNV clustering complete.\n');

clusterLabelsMsc = kmeans_simple(mscSpectra, nKmeansClusters, kmeansMaxIterations);
clusterMsc = reconstruct_label_image(finalMask, clusterLabelsMsc);
fprintf('  MSC clustering complete.\n');

% Visualize k-means comparison
figure('Name', 'K-means clustering comparison', 'Color', 'w');

subplot(2, 2, 1);
imshow(rgbImage);
title('Pseudo-RGB composite');

subplot(2, 2, 2);
imagesc(clusterRaw);
axis image off;
colorbar;
title('Raw clustering');

subplot(2, 2, 3);
imagesc(clusterSnv);
axis image off;
colorbar;
title('SNV clustering');

subplot(2, 2, 4);
imagesc(clusterMsc);
axis image off;
colorbar;
title('MSC clustering');

maybe_save_figure(gcf, 'clustering_comparison.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 15. Bean-level full-spectrum dataset for classification
%% =========================================================================
%% Extract mean full-range spectra for each bean from the ground-truth data.
%% These represent the aggregate spectral signature of each bean sample.

fprintf('\nExtracting bean-level full-spectrum dataset...\n');

fullSpectra = collect_matrix(beanData, 'Full_mean_spec');
fprintf('  Full-range spectra: %d x %d (beans x bands)\n', ...
    size(fullSpectra, 1), size(fullSpectra, 2));

% Visualize mean spectra by category
newMask = summaryTable.IsNew;
oldMask = ~summaryTable.IsNew;
germMask = summaryTable.Germination == 1;
nonGermMask = summaryTable.Germination == 0;

figure('Name', 'Mean spectra by category', 'Color', 'w');

subplot(2, 1, 1);
plot(fullWavelengths, mean(fullSpectra(newMask, :), 1), 'b-', 'LineWidth', 2, 'DisplayName', 'New');
hold on;
plot(fullWavelengths, mean(fullSpectra(oldMask, :), 1), 'r-', 'LineWidth', 2, 'DisplayName', 'Old');
hold off;
grid on;
xlabel('Wavelength (nm)');
ylabel('Reflectance');
legend;
title('Mean spectra by age');

subplot(2, 1, 2);
plot(fullWavelengths, mean(fullSpectra(germMask, :), 1), 'g-', 'LineWidth', 2, 'DisplayName', 'Germinated');
hold on;
plot(fullWavelengths, mean(fullSpectra(nonGermMask, :), 1), 'm-', 'LineWidth', 2, 'DisplayName', 'Non-germinated');
hold off;
grid on;
xlabel('Wavelength (nm)');
ylabel('Reflectance');
legend;
title('Mean spectra by germination');

maybe_save_figure(gcf, 'mean_spectra.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 16. Create classification labels and data split
%% =========================================================================
%% Two classification targets are defined:
%% 1. Age: new (1) vs old (0)
%% 2. Germination: germinated (1) vs non-germinated (0)
%%
%% A balanced calibration/test split is created. Note: the same split is
%% currently used for both age and germination classification. This is
%% acceptable only if the class balance remains valid for both targets.
%% If germination labels have a significantly different class distribution,
%% consider creating a separate balanced split for germination.

fprintf('\nPreparing classification dataset...\n');

ageLabels = double(summaryTable.IsNew);
germinationLabels = double(summaryTable.Germination);

fprintf('  Age labels: %d class-0 (old), %d class-1 (new)\n', ...
    sum(ageLabels == 0), sum(ageLabels == 1));
fprintf('  Germination labels: %d class-0 (non-germ), %d class-1 (germ)\n', ...
    sum(germinationLabels == 0), sum(germinationLabels == 1));

split = balanced_split(ageLabels, calibrationPerClass, testPerClass);
fprintf('  Calibration samples: %d (per class: %d)\n', ...
    numel(split.calibration), calibrationPerClass);
fprintf('  Test samples: %d (per class: %d)\n', ...
    numel(split.test), testPerClass);

fprintf('  NOTE: Same split used for both age and germination targets.\n');

%% =========================================================================
%% 17. PLS-DA age classification
%% =========================================================================
%% Partial Least Squares Discriminant Analysis (PLS-DA) is performed to
%% predict bean age (new vs old) from full-range spectra.
%% Three preprocessing variants are compared: raw, SNV, MSC.
%% Leave-one-out cross-validation is used to select the optimal number of
%% latent variables.

fprintf('\n========== AGE CLASSIFICATION (New vs Old) ==========\n');

agePlsda = struct();

for prep_idx = 1:3
    prep_names = {'raw', 'snv', 'msc'};
    prep_name = prep_names{prep_idx};
    
    fprintf('\nPreprocessing: %s\n', upper(prep_name));
    
    % Apply preprocessing
    Xprocessed = apply_preprocessing(fullSpectra, prep_name);
    
    % Split data
    Xcal = Xprocessed(split.calibration, :);
    ycal = ageLabels(split.calibration);
    Xtest = Xprocessed(split.test, :);
    ytest = ageLabels(split.test);
    
    % Leave-one-out cross-validation for component selection
    maxComponents = min([maxPlsComponentsCap, size(Xcal, 1) - 1, size(Xcal, 2)]);
    [bestComponents, cvAccuracy] = plsda_loocv(Xcal, ycal, maxComponents);
    
    fprintf('  Best components (LOOCV): %d, CV accuracy: %.4f\n', bestComponents, cvAccuracy);
    
    % Fit final model
    model = pls1_fit(Xcal, ycal, bestComponents);
    
    % Predict on test set
    yhat = pls1_predict(model, Xtest);
    predicted = double(yhat >= plsdaThreshold);
    conf = confusion_matrix(ytest, predicted);
    
    testAccuracy = sum(diag(conf)) / sum(conf(:));
    fprintf('  Test accuracy: %.4f\n', testAccuracy);
    fprintf('  Confusion matrix:\n');
    disp(conf);
    
    % Store results
    agePlsda.(prep_name) = struct(...
        'Name', 'Age', ...
        'Preprocessing', prep_name, ...
        'Components', bestComponents, ...
        'CVAccuracy', cvAccuracy, ...
        'TestScores', yhat, ...
        'PredictedLabels', predicted, ...
        'TrueLabels', ytest, ...
        'ConfusionMatrix', conf, ...
        'TestAccuracy', testAccuracy, ...
        'Wavelengths', fullWavelengths(:)');
end

fprintf('\n========== AGE CLASSIFICATION SUMMARY ==========\n');
fprintf('Raw:   CV acc: %.4f, Test acc: %.4f (%d LV)\n', ...
    agePlsda.raw.CVAccuracy, agePlsda.raw.TestAccuracy, agePlsda.raw.Components);
fprintf('SNV:   CV acc: %.4f, Test acc: %.4f (%d LV)\n', ...
    agePlsda.snv.CVAccuracy, agePlsda.snv.TestAccuracy, agePlsda.snv.Components);
fprintf('MSC:   CV acc: %.4f, Test acc: %.4f (%d LV)\n', ...
    agePlsda.msc.CVAccuracy, agePlsda.msc.TestAccuracy, agePlsda.msc.Components);

%% =========================================================================
%% 18. PLS-DA germination classification
%% =========================================================================
%% Partial Least Squares Discriminant Analysis (PLS-DA) is performed to
%% predict bean germination status (germinated vs non-germinated) from
%% full-range spectra, using the same methodology as the age classification.

fprintf('\n========== GERMINATION CLASSIFICATION (Germinated vs Non-germinated) ==========\n');

germinationPlsda = struct();

for prep_idx = 1:3
    prep_names = {'raw', 'snv', 'msc'};
    prep_name = prep_names{prep_idx};
    
    fprintf('\nPreprocessing: %s\n', upper(prep_name));
    
    % Apply preprocessing
    Xprocessed = apply_preprocessing(fullSpectra, prep_name);
    
    % Split data
    Xcal = Xprocessed(split.calibration, :);
    ycal = germinationLabels(split.calibration);
    Xtest = Xprocessed(split.test, :);
    ytest = germinationLabels(split.test);
    
    % Leave-one-out cross-validation for component selection
    maxComponents = min([maxPlsComponentsCap, size(Xcal, 1) - 1, size(Xcal, 2)]);
    [bestComponents, cvAccuracy] = plsda_loocv(Xcal, ycal, maxComponents);
    
    fprintf('  Best components (LOOCV): %d, CV accuracy: %.4f\n', bestComponents, cvAccuracy);
    
    % Fit final model
    model = pls1_fit(Xcal, ycal, bestComponents);
    
    % Predict on test set
    yhat = pls1_predict(model, Xtest);
    predicted = double(yhat >= plsdaThreshold);
    conf = confusion_matrix(ytest, predicted);
    
    testAccuracy = sum(diag(conf)) / sum(conf(:));
    fprintf('  Test accuracy: %.4f\n', testAccuracy);
    fprintf('  Confusion matrix:\n');
    disp(conf);
    
    % Store results
    germinationPlsda.(prep_name) = struct(...
        'Name', 'Germination', ...
        'Preprocessing', prep_name, ...
        'Components', bestComponents, ...
        'CVAccuracy', cvAccuracy, ...
        'TestScores', yhat, ...
        'PredictedLabels', predicted, ...
        'TrueLabels', ytest, ...
        'ConfusionMatrix', conf, ...
        'TestAccuracy', testAccuracy, ...
        'Wavelengths', fullWavelengths(:)');
end

fprintf('\n========== GERMINATION CLASSIFICATION SUMMARY ==========\n');
fprintf('Raw:   CV acc: %.4f, Test acc: %.4f (%d LV)\n', ...
    germinationPlsda.raw.CVAccuracy, germinationPlsda.raw.TestAccuracy, germinationPlsda.raw.Components);
fprintf('SNV:   CV acc: %.4f, Test acc: %.4f (%d LV)\n', ...
    germinationPlsda.snv.CVAccuracy, germinationPlsda.snv.TestAccuracy, germinationPlsda.snv.Components);
fprintf('MSC:   CV acc: %.4f, Test acc: %.4f (%d LV)\n', ...
    germinationPlsda.msc.CVAccuracy, germinationPlsda.msc.TestAccuracy, germinationPlsda.msc.Components);

%% =========================================================================
%% 19. Summary figures
%% =========================================================================
%% Generate a comprehensive figure showing PLS-DA confusion matrices
%% for both classification targets across all preprocessing methods.

fprintf('\nGenerating PLS-DA summary figure...\n');

figure('Name', 'PLS-DA summary', 'Color', 'w');

% Age classification
for prep_idx = 1:3
    prep_names = {'raw', 'snv', 'msc'};
    prep_name = prep_names{prep_idx};
    item = agePlsda.(prep_name);
    
    subplot(2, 3, prep_idx);
    imagesc(item.ConfusionMatrix);
    axis image;
    colorbar;
    title(sprintf('Age (%s)\nAcc: %.2f | LV: %d', upper(prep_name), item.TestAccuracy, item.Components));
    xlabel('Predicted');
    ylabel('True');
    caxis([0, max(item.ConfusionMatrix(:))]);
end

% Germination classification
for prep_idx = 1:3
    prep_names = {'raw', 'snv', 'msc'};
    prep_name = prep_names{prep_idx};
    item = germinationPlsda.(prep_name);
    
    subplot(2, 3, 3 + prep_idx);
    imagesc(item.ConfusionMatrix);
    axis image;
    colorbar;
    title(sprintf('Germination (%s)\nAcc: %.2f | LV: %d', upper(prep_name), item.TestAccuracy, item.Components));
    xlabel('Predicted');
    ylabel('True');
    caxis([0, max(item.ConfusionMatrix(:))]);
end

maybe_save_figure(gcf, 'plsda_summary.png', saveFigures, closeFiguresAfterSaving, outputDir);

%% =========================================================================
%% 20. Results structure and output saving
%% =========================================================================
%% Assemble all major analysis artifacts into a single `results` structure
%% and optionally save to MAT file.

fprintf('\nAssembling results structure...\n');

results = struct();

% Metadata
results.summaryTable = summaryTable;
results.vnirWavelengths = vnirWavelengths;
results.fullWavelengths = fullWavelengths;

% Initial processing
results.rgbImage = rgbImage;

% Masks and thresholds
results.ratioImage = ratioImage;
results.ratioThreshold = ratioThreshold;
results.firstMask = firstMask;
results.pc2ScoreImage = pc2ScoreImage;
results.pc2Threshold = pc2Threshold;
results.secondMask = secondMask;
results.pc3ScoreImage = pc3ScoreImage;
results.pc3Threshold = pc3Threshold;
results.thirdMask = thirdMask;
results.finalMask = finalMask;

% First-pass PCA
results.firstPassPCA = struct(...
    'Loadings', loadingsFirst, ...
    'Scores', scoresFirst, ...
    'Eigenvalues', eigenvaluesFirst, ...
    'ScoreImages', scoreImagesFirst);

% Spectra
results.maskedSpectra = struct(...
    'raw', maskedSpectra, ...
    'snv', snvSpectra, ...
    'msc', mscSpectra);

% PCA for all variants
results.pca = struct(...
    'raw', struct('Loadings', loadingsRaw, 'Scores', scoresRaw, 'Eigenvalues', eigenvaluesRaw, 'ScoreImages', scoreImagesRaw), ...
    'snv', struct('Loadings', loadingsSnv, 'Scores', scoresSnv, 'Eigenvalues', eigenvaluesSnv, 'ScoreImages', scoreImagesSnv), ...
    'msc', struct('Loadings', loadingsMsc, 'Scores', scoresMsc, 'Eigenvalues', eigenvaluesMsc, 'ScoreImages', scoreImagesMsc));

% Clustering
results.clustering = struct(...
    'raw', clusterRaw, ...
    'snv', clusterSnv, ...
    'msc', clusterMsc);

% Classification
results.classification = struct(...
    'fullSpectra', fullSpectra, ...
    'ageLabels', ageLabels, ...
    'germinationLabels', germinationLabels, ...
    'split', split, ...
    'agePlsda', agePlsda, ...
    'germinationPlsda', germinationPlsda);

% Parameters used
results.parameters = struct(...
    'rgbWavelengths', rgbWavelengths, ...
    'ratioHighWavelength', ratioHighWavelength, ...
    'ratioLowWavelength', ratioLowWavelength, ...
    'nFirstPassPCs', nFirstPassPCs, ...
    'pcForSecondMask', pcForSecondMask, ...
    'pcForThirdMask', pcForThirdMask, ...
    'nKmeansClusters', nKmeansClusters, ...
    'calibrationPerClass', calibrationPerClass, ...
    'testPerClass', testPerClass, ...
    'plsdaThreshold', plsdaThreshold);

fprintf('Results structure assembled.\n');

if saveResults
    resultsFilePath = fullfile(outputDir, 'anon_thesis_workflow_results.mat');
    save(resultsFilePath, 'results', '-v7.3');
    fprintf('\nResults saved to: %s\n', resultsFilePath);
end

fprintf('\n========== WORKFLOW COMPLETE ==========\n\n');

%% =========================================================================
%% LOCAL HELPER FUNCTIONS
%% =========================================================================

%% Build summary table from bean data
function summaryTable = build_summary_table(beanData)
numBeans = numel(beanData);
names = strings(numBeans, 1);
isNew = false(numBeans, 1);
germination = zeros(numBeans, 1);

for idx = 1:numBeans
    names(idx) = string(beanData(idx).Name);
    isNew(idx) = contains(lower(names(idx)), 'new');
    germination(idx) = double(beanData(idx).Germination);
end

summaryTable = table((1:numBeans)', names, isNew, germination, ...
    'VariableNames', {'BeanIndex', 'Name', 'IsNew', 'Germination'});
end

%% Find band index closest to target wavelength
function index = nearest_band(wavelengths, target)
[~, index] = min(abs(double(wavelengths(:)) - target));
end

%% Scale array to [0, 1]
function scaled = scale_01(X)
X = double(X);
minValue = min(X(:));
maxValue = max(X(:));
if maxValue <= minValue
    scaled = zeros(size(X));
else
    scaled = (X - minValue) ./ (maxValue - minValue);
end
end

%% Otsu automatic threshold
function threshold = otsu_threshold(values)
values = double(values(:));
values = values(~isnan(values) & ~isinf(values));
values = scale_01(values);

if isempty(values)
    threshold = 0.5;
    return;
end

numBins = 256;
edges = linspace(0, 1, numBins + 1);
counts = histcounts(values, edges);
probabilities = counts / sum(counts);
binCenters = (edges(1:end-1) + edges(2:end)) / 2;

omega = cumsum(probabilities);
mu = cumsum(probabilities .* binCenters);
muT = mu(end);

sigmaBetween = (muT * omega - mu) .^ 2 ./ max(omega .* (1 - omega), eps);
sigmaBetween(~isfinite(sigmaBetween)) = 0;
[~, maxIndex] = max(sigmaBetween);
threshold = binCenters(maxIndex);
end

%% Extract spectra from masked region
function maskedSpectra = extract_masked_spectra(cube, mask)
[rows, cols, bands] = size(cube);
unfoldedCube = reshape(cube, rows * cols, bands);
maskedSpectra = unfoldedCube(mask(:), :);
end

%% Reconstruct score images from mask and scores
function scoreImages = reconstruct_score_images(mask, scores, nComponents)
scoreImages = zeros(size(mask, 1), size(mask, 2), nComponents);
for componentIdx = 1:nComponents
    plane = nan(size(mask));
    plane(mask) = scores(:, componentIdx);
    scoreImages(:, :, componentIdx) = plane;
end
end

%% Create score-based mask using Otsu threshold
function [refinedMask, threshold] = score_mask(scoreImage, supportMask, centralRegionFraction, ...
    minimumScorePixels, minimumRetainedPixelsAbsolute, minimumRetainedPixelsFraction)

validScores = scoreImage(supportMask);
validScores = validScores(~isnan(validScores));

if numel(validScores) < minimumScorePixels
    refinedMask = supportMask;
    threshold = nan;
    return;
end

threshold = otsu_threshold(validScores);
maskA = supportMask & (scoreImage >= threshold);
maskB = supportMask & (scoreImage < threshold);

% Choose polarity based on central support
centerRows = round(size(scoreImage, 1) * centralRegionFraction(1)) : ...
             round(size(scoreImage, 1) * centralRegionFraction(2));
centerCols = round(size(scoreImage, 2) * centralRegionFraction(1)) : ...
             round(size(scoreImage, 2) * centralRegionFraction(2));
centerSupport = false(size(scoreImage));
centerSupport(centerRows, centerCols) = true;

if nnz(maskA & centerSupport) >= nnz(maskB & centerSupport)
    refinedMask = maskA;
else
    refinedMask = maskB;
end

% Fallback if refinement is too aggressive
minRetainedPixels = max(minimumRetainedPixelsAbsolute, ...
    round(minimumRetainedPixelsFraction * nnz(supportMask)));
if nnz(refinedMask) < minRetainedPixels
    refinedMask = supportMask;
end
end

%% MSC preprocessing
function corrected = msc(X)
reference = mean(X, 1);
corrected = zeros(size(X));
refAugmented = [reference(:), ones(numel(reference), 1)];

for rowIdx = 1:size(X, 1)
    coeffs = refAugmented \ X(rowIdx, :)';
    slope = coeffs(1);
    intercept = coeffs(2);
    corrected(rowIdx, :) = (X(rowIdx, :) - intercept) ./ slope;
end
end

%% Simple k-means clustering
function labels = kmeans_simple(X, nClusters, maxIterations)
rng(7);
numSamples = size(X, 1);
seedIdx = randperm(numSamples, nClusters);
centers = X(seedIdx, :);
labels = ones(numSamples, 1);

for iterIdx = 1:maxIterations
    distances = zeros(numSamples, nClusters);
    for clusterIdx = 1:nClusters
        deltas = X - centers(clusterIdx, :);
        distances(:, clusterIdx) = sum(deltas .^ 2, 2);
    end
    
    [~, newLabels] = min(distances, [], 2);
    if isequal(newLabels, labels)
        break;
    end
    labels = newLabels;
    
    for clusterIdx = 1:nClusters
        members = X(labels == clusterIdx, :);
        if isempty(members)
            centers(clusterIdx, :) = X(randi(numSamples), :);
        else
            centers(clusterIdx, :) = mean(members, 1);
        end
    end
end
end

%% Reconstruct label image from mask and labels
function labelImage = reconstruct_label_image(mask, labels)
labelImage = nan(size(mask));
labelImage(mask) = labels;
end

%% Collect a field from bean data structure into a matrix
function matrix = collect_matrix(beanData, fieldName)
numBeans = numel(beanData);
template = beanData(1).(fieldName);
matrix = zeros(numBeans, numel(template));
for idx = 1:numBeans
    matrix(idx, :) = reshape(double(beanData(idx).(fieldName)), 1, []);
end
end

%% Create balanced calibration/test split
function split = balanced_split(labels, calibrationPerClass, testPerClass)
classOne = find(labels == 1);
classZero = find(labels == 0);

minimumClassCount = min(numel(classOne), numel(classZero));
if minimumClassCount <= calibrationPerClass
    error('Not enough samples to create the requested balanced split.');
end

effectiveTestPerClass = min(testPerClass, minimumClassCount - calibrationPerClass);
if effectiveTestPerClass < testPerClass
    warning('Requested %d test samples per class, using %d based on available data.', ...
        testPerClass, effectiveTestPerClass);
end

split.calibration = [classOne(1:calibrationPerClass); classZero(1:calibrationPerClass)];
split.test = [classOne(calibrationPerClass + (1:effectiveTestPerClass)); ...
    classZero(calibrationPerClass + (1:effectiveTestPerClass))];
end

%% Apply preprocessing to spectral data
function Xprocessed = apply_preprocessing(X, preprocessingName)
switch lower(preprocessingName)
    case 'raw'
        Xprocessed = X;
    case 'snv'
        Xprocessed = SNV(X);
    case 'msc'
        Xprocessed = msc(X);
    otherwise
        error('Unknown preprocessing option: %s', preprocessingName);
end
end

%% Leave-one-out cross-validation for PLS-DA component selection
function [bestComponents, bestAccuracy] = plsda_loocv(X, y, maxComponents)
numSamples = size(X, 1);
accuracies = zeros(maxComponents, 1);

for componentIdx = 1:maxComponents
    predictions = zeros(numSamples, 1);
    for sampleIdx = 1:numSamples
        calibrationMask = true(numSamples, 1);
        calibrationMask(sampleIdx) = false;
        model = pls1_fit(X(calibrationMask, :), y(calibrationMask), componentIdx);
        predictions(sampleIdx) = pls1_predict(model, X(sampleIdx, :));
    end
    predictedLabels = double(predictions >= 0.5);
    accuracies(componentIdx) = mean(predictedLabels == y);
end

[bestAccuracy, bestComponents] = max(accuracies);
end

%% PLS1 model fitting
function model = pls1_fit(X, y, nComponents)
X = double(X);
y = double(y(:));

xMean = mean(X, 1);
yMean = mean(y);
E = X - xMean;
f = y - yMean;

[numSamples, numVariables] = size(E);
W = zeros(numVariables, nComponents);
P = zeros(numVariables, nComponents);
Q = zeros(nComponents, 1);
T = zeros(numSamples, nComponents);

for componentIdx = 1:nComponents
    w = E' * f;
    wNorm = norm(w);
    if wNorm == 0
        break;
    end
    w = w / wNorm;
    t = E * w;
    denom = t' * t;
    if denom == 0
        break;
    end
    p = (E' * t) / denom;
    q = (f' * t) / denom;
    
    E = E - t * p';
    f = f - t * q;
    
    W(:, componentIdx) = w;
    P(:, componentIdx) = p;
    Q(componentIdx) = q;
    T(:, componentIdx) = t;
end

validComponents = find(any(W ~= 0, 1));
W = W(:, validComponents);
P = P(:, validComponents);
Q = Q(validComponents);
T = T(:, validComponents);

beta = W / (P' * W) * Q;
intercept = yMean - xMean * beta;

model = struct('Beta', beta, 'Intercept', intercept, 'XMean', xMean, 'YMean', yMean, 'Scores', T);
end

%% PLS1 model prediction
function yhat = pls1_predict(model, X)
X = double(X);
yhat = X * model.Beta + model.Intercept;
end

%% Confusion matrix
function confusion = confusion_matrix(yTrue, yPred)
confusion = zeros(2, 2);
for idx = 1:numel(yTrue)
    confusion(yTrue(idx) + 1, yPred(idx) + 1) = confusion(yTrue(idx) + 1, yPred(idx) + 1) + 1;
end
end

%% Conditionally save figure
function maybe_save_figure(figHandle, fileName, saveFigures, closeFiguresAfterSaving, outputDir)
if saveFigures
    if isempty(outputDir)
        error('Cannot save figures: outputDir is empty.');
    end
    filePath = fullfile(outputDir, fileName);
    saveas(figHandle, filePath);
end

if closeFiguresAfterSaving
    close(figHandle);
end
end
