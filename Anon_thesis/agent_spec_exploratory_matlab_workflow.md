# Coding Agent Specification: Rework `run_anon_thesis_workflow.m` into an Exploratory MATLAB Script

## 1. Objective

Rework the current MATLAB file:

```matlab
run_anon_thesis_workflow.m
```

into a new exploratory analysis script:

```matlab
explore_anon_thesis_workflow.m
```

The goal is **not** to redesign the scientific method. The goal is to preserve the current analytical logic while changing the structure from a function-heavy automated pipeline into a transparent, section-by-section exploratory MATLAB script.

The reworked script should be suitable for thesis-stage analysis, where the user needs to run sections manually, inspect figures, adjust thresholds, inspect intermediate variables, and continue.

The current script is too automated for exploratory development: it is wrapped as a function, hides many decisions inside local helper functions, saves figures automatically, closes figures immediately, and returns a final `results` structure. The new script should keep the same scientific flow but expose the workflow more clearly.

---

## 2. High-level Instruction

Convert the current workflow from this style:

```matlab
function results = run_anon_thesis_workflow()
    % automated pipeline
    % many local helper functions
    % saves outputs
    % closes figures
end
```

into this style:

```matlab
%% 1. Setup and file paths
%% 2. Load data
%% 3. Initial image inspection
%% 4. First-pass mask
%% 5. First-pass PCA
%% 6. PCA-guided mask refinement
%% 7. Final masked spectra extraction
%% 8. Preprocessing comparison
%% 9. PCA comparison
%% 10. K-means clustering
%% 11. Bean-level spectra extraction
%% 12. Calibration/test split
%% 13. PLS-DA age classification
%% 14. PLS-DA germination classification
%% 15. Summary figures
%% 16. Results structure and saving
```

Important intermediate variables should remain visible in the MATLAB workspace.

Do **not** convert the workflow into a fully modular production pipeline. The user wants an exploratory script, not a black-box function.

---

## 3. Preserve the Existing Scientific Workflow

The new script must preserve the current analytical sequence.

The current script performs the following operations:

1. Defines base and output directories.
2. Loads the VNIR cube from `nb_s6_vnir_2.mat`.
3. Checks that the loaded variable is named `im`.
4. Converts the image cube to `double`.
5. Loads bean-level ground-truth data from `kidney_bean_data_GT.mat`.
6. Checks that the loaded variable is named `kidney_bean_data`.
7. Converts the bean data structure into a row vector.
8. Builds a summary table with:
   - bean index
   - bean name
   - new/old label
   - germination label
9. Extracts VNIR wavelengths from `beanData(1).VNIR_WL`.
10. Extracts full-range wavelengths from `beanData(1).Full_WL`.
11. Creates a pseudo-RGB image using bands nearest to:
   - 650 nm
   - 550 nm
   - 450 nm
12. Creates a first-pass ratio image using:
   - 900 nm band
   - 650 nm band
13. Scales the ratio image to 0–1.
14. Applies an Otsu threshold to create a first-pass mask.
15. Chooses the mask polarity based on support in the central image region.
16. Extracts masked spectra from the first-pass mask.
17. Runs PCA on the first-pass masked spectra.
18. Reconstructs the first three PCA score images.
19. Uses PC2 to create a second mask.
20. Uses PC3 to create a third mask.
21. Combines the first, second, and third masks into a final mask.
22. Extracts final masked spectra.
23. Applies SNV preprocessing.
24. Applies MSC preprocessing.
25. Runs PCA separately on:
   - raw masked spectra
   - SNV spectra
   - MSC spectra
26. Reconstructs PCA score images for raw, SNV, and MSC data.
27. Runs k-means clustering with two clusters for:
   - raw spectra
   - SNV spectra
   - MSC spectra
28. Reconstructs cluster maps from masked-pixel labels.
29. Extracts bean-level full spectra from `Full_mean_spec`.
30. Creates age labels from `IsNew`.
31. Creates germination labels from `Germination`.
32. Creates a balanced calibration/test split using:
   - 75 calibration samples per class
   - 25 test samples per class
33. Runs PLS-DA classification for age.
34. Runs PLS-DA classification for germination.
35. Compares raw, SNV, and MSC preprocessing for both classification tasks.
36. Uses leave-one-out cross-validation to choose the number of latent variables.
37. Uses a classification threshold of `0.5`.
38. Generates figures for:
   - RGB and mask workflow
   - PCA bundles
   - k-means comparison
   - mean spectra
   - PLS-DA summary
39. Builds a final `results` structure.
40. Saves the results to:

```matlab
analysis_outputs/anon_thesis_workflow_results.mat
```

The reworked script must keep this workflow unless there is a clear coding error.

---

## 4. Target Script Structure

Use MATLAB section breaks. The script should be readable when run one section at a time.

Suggested structure:

```matlab
%% 1. Setup and file paths

%% 2. Load VNIR image cube

%% 3. Load bean-level ground-truth data

%% 4. Build bean summary table

%% 5. Initial RGB / pseudo-RGB visualisation

%% 6. First-pass 900/650 nm ratio mask

%% 7. First-pass PCA on masked spectra

%% 8. PCA-guided mask refinement using PC2 and PC3

%% 9. Final masked spectra extraction

%% 10. Spectral preprocessing: raw, SNV, MSC

%% 11. PCA comparison across preprocessing methods

%% 12. K-means clustering comparison

%% 13. Bean-level full-spectrum dataset for classification

%% 14. Balanced calibration/test split

%% 15. PLS-DA age classification

%% 16. PLS-DA germination classification

%% 17. Summary figures

%% 18. Results structure and output saving

%% Local helper functions
```

The script can still include helper functions at the bottom, but the main workflow should be readable from top to bottom without needing to jump constantly between functions.

---

## 5. Key Exploratory Parameters to Expose

Move important choices out of hidden helper functions and place them in visible parameter blocks near the top of the script or near the relevant section.

At minimum, expose these parameters:

```matlab
rgbWavelengths = [650, 550, 450];

ratioHighWavelength = 900;
ratioLowWavelength = 650;

nFirstPassPCs = 3;
pcForSecondMask = 2;
pcForThirdMask = 3;

centralRegionFraction = [0.3, 0.7];

minimumScorePixels = 10;
minimumRetainedPixelsAbsolute = 50;
minimumRetainedPixelsFraction = 0.2;

nKmeansClusters = 2;
kmeansMaxIterations = 50;
randomSeed = 7;

calibrationPerClass = 75;
testPerClass = 25;

maxPlsComponentsCap = 20;
plsdaThreshold = 0.5;

saveFigures = true;
closeFiguresAfterSaving = false;
saveResults = true;
```

Thresholds calculated by Otsu should also be stored in visible variables, for example:

```matlab
ratioThreshold
pc2Threshold
pc3Threshold
```

If helper functions calculate thresholds, the thresholds should be returned and stored.

---

## 6. Important Structural Changes

### 6.1 Remove the Function Wrapper

The new file should be a script, not a function.

Do not start with:

```matlab
function results = run_anon_thesis_workflow()
```

Instead, start with a normal script setup section:

```matlab
%% 1. Setup and file paths

clearvars;
close all;
clc;
```

If the user may want to keep existing workspace variables, make `clearvars` optional:

```matlab
resetWorkspace = true;

if resetWorkspace
    clearvars;
    close all;
    clc;
end
```

Because this is an exploratory script, avoid being overly aggressive with clearing variables unless clearly documented.

---

### 6.2 Keep Intermediate Variables Visible

The following variables should exist in the workspace after relevant sections are run:

```matlab
cube
beanData
summaryTable
vnirWavelengths
fullWavelengths

rgbImage
ratioImage
ratioThreshold
firstMask

maskedSpectraFirst
loadingsFirst
scoresFirst
eigenvaluesFirst
scoreImagesFirst

pc2ScoreImage
pc3ScoreImage
pc2Threshold
pc3Threshold
secondMask
thirdMask
finalMask

maskedSpectra
snvSpectra
mscSpectra

loadingsRaw
scoresRaw
eigenvaluesRaw
scoreImagesRaw

loadingsSnv
scoresSnv
eigenvaluesSnv
scoreImagesSnv

loadingsMsc
scoresMsc
eigenvaluesMsc
scoreImagesMsc

clusterLabelsRaw
clusterLabelsSnv
clusterLabelsMsc
clusterRaw
clusterSnv
clusterMsc

fullSpectra
ageLabels
germinationLabels
split

agePlsda
germinationPlsda
results
```

Do not hide these inside a `results` structure only.

---

### 6.3 Expose Masking Logic

The current masking logic is methodologically important and should be visible in the main script.

The script should show the ratio mask logic directly:

```matlab
highIdx = nearest_band(vnirWavelengths, ratioHighWavelength);
lowIdx = nearest_band(vnirWavelengths, ratioLowWavelength);

ratioImage = cube(:, :, highIdx) ./ max(cube(:, :, lowIdx), eps);
ratioImage = scale_01(ratioImage);

ratioThreshold = otsu_threshold(ratioImage(:));

maskA = ratioImage >= ratioThreshold;
maskB = ratioImage < ratioThreshold;

% Choose polarity based on central support.
```

The script should also show PC-based mask refinement explicitly:

```matlab
pc2ScoreImage = scoreImagesFirst(:, :, pcForSecondMask);
[secondMask, pc2Threshold] = score_mask(pc2ScoreImage, firstMask, ...);

pc3ScoreImage = scoreImagesFirst(:, :, pcForThirdMask);
[thirdMask, pc3Threshold] = score_mask(pc3ScoreImage, firstMask & secondMask, ...);

finalMask = firstMask & secondMask & thirdMask;
```

Store and plot all masks.

---

### 6.4 Keep Diagnostic Figures Open by Default

The current script saves figures and immediately closes them:

```matlab
saveas(gcf, fullfile(outputDir, 'masking_workflow.png'));
close(gcf);
```

For exploratory use, figures should remain open by default.

Use this pattern:

```matlab
if saveFigures
    saveas(gcf, fullfile(outputDir, 'masking_workflow.png'));
end

if closeFiguresAfterSaving
    close(gcf);
end
```

Set:

```matlab
closeFiguresAfterSaving = false;
```

by default.

---

### 6.5 Keep Helper Functions, but Only for Low-level Operations

Acceptable helper functions at the bottom of the script:

```matlab
nearest_band
scale_01
otsu_threshold
extract_masked_spectra
reconstruct_score_images
reconstruct_label_image
msc
apply_preprocessing
kmeans_simple
balanced_split
plsda_loocv
pls1_fit
pls1_predict
confusion_matrix
maybe_save_figure
```

Avoid hiding the scientific flow inside large helper functions such as:

```matlab
local_run_plsda_suite
local_plot_pca_bundle
local_plot_plsda_bundle
```

These can be kept only if the main script clearly shows the inputs, outputs, and interpretation. Prefer expanding the workflow in the main script for readability.

---

## 7. Visualisation Requirements

The script should generate diagnostic figures after major transformations.

### 7.1 Initial Image Inspection

Show:

```matlab
imshow(rgbImage)
title('Pseudo-RGB composite')
```

Also consider showing individual bands near 450, 550, 650, and 900 nm.

---

### 7.2 Ratio Mask Figure

Show:

1. RGB image
2. Ratio image
3. Ratio histogram with Otsu threshold
4. First-pass mask

---

### 7.3 First-pass PCA Figure

Show:

1. Explained variance for first-pass PCA
2. PC1 score image
3. PC2 score image
4. PC3 score image
5. PC1 loading
6. PC2 loading
7. PC3 loading

---

### 7.4 Mask Refinement Figure

Show:

1. First-pass mask
2. PC2 score image
3. Second mask
4. PC3 score image
5. Third mask
6. Final mask

---

### 7.5 Preprocessing Comparison Figure

Show mean spectra or example spectra for:

1. Raw masked spectra
2. SNV spectra
3. MSC spectra

---

### 7.6 PCA Comparison Figure

For raw, SNV, and MSC data, show:

1. Explained variance
2. PC1 loading
3. PC2 loading
4. PC3 loading
5. PC1 score image
6. PC2 score image

This can be done in separate figures for each preprocessing method.

---

### 7.7 K-means Figure

Show:

1. RGB image
2. Raw clustering map
3. SNV clustering map
4. MSC clustering map

---

### 7.8 Mean Spectra Figure

Show:

1. Mean full-range spectra for new vs old beans
2. Mean full-range spectra for germinated vs non-germinated beans

---

### 7.9 PLS-DA Summary Figure

Show confusion matrices for:

1. Age classification, raw
2. Age classification, SNV
3. Age classification, MSC
4. Germination classification, raw
5. Germination classification, SNV
6. Germination classification, MSC

---

## 8. Classification Workflow Requirements

Preserve the current PLS-DA logic.

The script should still perform PLS-DA for:

```matlab
ageLabels
germinationLabels
```

using:

```matlab
fullSpectra = collect_matrix(beanData, 'Full_mean_spec');
```

Preserve comparison across:

```matlab
raw
snv
msc
```

Preserve leave-one-out cross-validation for selecting the best number of latent variables.

Preserve the threshold:

```matlab
predictedLabels = double(yhat >= plsdaThreshold);
```

where:

```matlab
plsdaThreshold = 0.5;
```

Make the chosen number of latent variables, cross-validation accuracy, test scores, predicted labels, true labels, confusion matrix, and test accuracy visible.

For each classification target and preprocessing method, keep structures such as:

```matlab
agePlsda.raw
agePlsda.snv
agePlsda.msc

germinationPlsda.raw
germinationPlsda.snv
germinationPlsda.msc
```

Each should contain:

```matlab
Name
Preprocessing
Components
CVAccuracy
TestScores
PredictedLabels
TrueLabels
ConfusionMatrix
TestAccuracy
Wavelengths
```

---

## 9. Data-splitting Requirements

Preserve the current balanced split approach:

```matlab
split = balanced_split(ageLabels, calibrationPerClass, testPerClass);
```

However, make an important note in the script comments:

```matlab
% The same split is currently used for age and germination classification.
% This is acceptable only if the class balance remains valid for both targets.
% If germination labels have a different class distribution, consider creating
% a separate balanced split for germination.
```

Do not silently change the split strategy unless instructed.

---

## 10. Naming Conventions

Use clearer names than the current local-function style.

Prefer:

```matlab
cube
beanData
summaryTable
vnirWavelengths
fullWavelengths

rgbImage
ratioImage
firstMask
secondMask
thirdMask
finalMask

maskedSpectraFirst
maskedSpectra
snvSpectra
mscSpectra

scoreImagesFirst
scoreImagesRaw
scoreImagesSnv
scoreImagesMsc

clusterRaw
clusterSnv
clusterMsc

fullSpectra
ageLabels
germinationLabels
agePlsda
germinationPlsda
```

Avoid excessive use of generic names like:

```matlab
X
y
idx
item
suite
```

inside the main script. These are acceptable inside small helper functions.

---

## 11. Coding Style

Use clear section headings and short explanatory comments.

Prefer comments that explain why a step exists, not only what MATLAB syntax does.

Good:

```matlab
% The 900/650 nm ratio is used as a first-pass contrast image because it
% separates bean tissue from background more clearly than a single band.
```

Less useful:

```matlab
% Divide high band by low band.
```

Avoid very long comment blocks inside simple algorithms unless needed. In the current file, the k-means function contains a very long explanatory block. Shorten this to a concise comment, because the exploratory script should be readable without excessive internal documentation.

---

## 12. Non-goals

Do **not** do the following unless explicitly requested:

1. Do not replace PLS-DA with LDA, SVM, random forest, neural networks, or another classifier.
2. Do not replace the masking strategy with a different segmentation method.
3. Do not remove PCA.
4. Do not remove SNV or MSC.
5. Do not remove k-means clustering.
6. Do not remove age classification.
7. Do not remove germination classification.
8. Do not turn the script into a set of many separate files.
9. Do not hide the workflow inside high-level functions.
10. Do not optimise for compactness at the expense of interpretability.
11. Do not automatically close all figures by default.
12. Do not remove final saving of the results structure.

---

## 13. Outputs

The reworked script should create or preserve:

```matlab
analysis_outputs/
```

and save:

```matlab
analysis_outputs/anon_thesis_workflow_results.mat
```

If `saveFigures = true`, save figures such as:

```matlab
masking_workflow.png
first_pass_pca.png
mask_refinement.png
pca_raw.png
pca_snv.png
pca_msc.png
kmeans_comparison.png
mean_spectra.png
plsda_summary.png
```

The script should still be useful even if `saveFigures = false`; in that case, figures should appear interactively but not be saved.

---

## 14. Suggested Acceptance Criteria

The reworked script is successful if:

1. It is a script, not a function.
2. It can be run section by section in MATLAB.
3. Major intermediate variables remain visible in the workspace.
4. The original analytical sequence is preserved.
5. The ratio mask, PC2 mask, PC3 mask, and final mask are all visible and saved in variables.
6. Raw, SNV, and MSC preprocessing are all retained.
7. PCA is performed for raw, SNV, and MSC spectra.
8. K-means maps are generated for raw, SNV, and MSC spectra.
9. PLS-DA is performed for both age and germination classification.
10. Confusion matrices and test accuracies are reported.
11. Figures remain open by default.
12. Results are saved at the end if `saveResults = true`.
13. The code does not rely on hidden workspace variables.
14. Important parameters can be adjusted from visible parameter blocks.
15. The script feels like a guided exploratory thesis workflow rather than an automated black-box pipeline.

---

## 15. Final Instruction to the Coding Agent

Rework `run_anon_thesis_workflow.m` into `explore_anon_thesis_workflow.m`.

Preserve the scientific logic, but expose the workflow. Keep the analysis exploratory, visual, and section-based. Make important intermediate variables visible. Keep figures open by default. Keep final saving of figures and results optional through parameters. Use helper functions only for low-level repeated operations, not to hide the main scientific workflow.

The final script should allow the user to inspect each stage of the analysis: data loading, RGB visualisation, ratio masking, PCA-guided mask refinement, preprocessing comparison, PCA comparison, k-means clustering, bean-level PLS-DA classification, and final result saving.
