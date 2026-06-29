clear vars
close all 
cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Analysis')



% Loading in the structure so we can work with the data 
S = load('Merged_VNIR_SWIR_Summary_With_Combined.mat');
fieldnames(S);
T = S.T_merged;

load('Tbl_vnir_mean_mix.mat', 'Tbl_vnir_mean_mix');
load('Tbl_swir_mean_mix.mat', 'Tbl_swir_mean_mix');
load('Tbl_combined_mix.mat', 'Tbl_combined_mix');



%% Plotting the mean spectra per class for the combined dataset (training)

WL_col = startsWith(Tbl_combined_mix.Properties.VariableNames, "WL_");
WL_names = Tbl_combined_mix.Properties.VariableNames(WL_col);
Wavelengths = str2double(extractAfter(WL_names, 'WL_'));

% Defining the x for which the mean will be calcualted 
X = Tbl_combined_mix{:, WL_col};
label = Tbl_combined_mix.Label; 

classes = categories(label); % defining the classes

figure;
hold on;

for i = 1:numel(classes)
    classNames = classes{i};
    classIdx = label == classNames;
    meanSpectrum = mean(X(classIdx, :), 1); 
    plot(Wavelengths, meanSpectrum, 'DisplayName', classNames);
end

xlabel('Wavelength (nm)');
ylabel('Mean Reflectance');
title('Mean Spectra per Class (Independent Test Set)');
legend('Location', 'best');
grid on;
hold off;

% reworking the table in real time 
outputFolder = 'C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Figures';

exportgraphics(gcf, fullfile(outputFolder, 'Mean_Spectra_Independent_Test_Set.png'), 'Resolution', 300);

% Doing the part with the PCA

% analysing the outliers of the dataset 




%%

% finding WL vector
vnirFile = 'BEN_VNIR.mat';
swirFile = 'PEN_SWIR.mat';


Sv = load(vnirFile);
Ss = load(swirFile);

WL_vnir = Sv.WL(:);
WL_swir = Ss.WL(:);

WL_combined = [WL_vnir; WL_swir(9:end)];

nVNIR = numel(WL_vnir);
nSWIR = numel(WL_swir) - 8;

clear Sv Ss

%% Preparing the dataset

% Converting cell spectra to numeric matrix 
specLen = cellfun(@numel, T.Combined_Spectrum);
assert(all(specLen == specLen(1)), "Combined_Spectrum lengths are inconsistent.");
X_raw = cell2mat(T.Combined_Spectrum);

X_vnir = X_raw(:, 1:nVNIR);
X_swir = X_raw(:, nVNIR+1 : nVNIR+nSWIR);

labels = extractBefore(string(T.Sample_ID), "_Sample_");
Y = categorical(labels);

% Keeping only valid rows
valid = all(isfinite(X_raw),2) & ~ismissing(Y);
X_raw = X_raw(valid,:);
X_vnir = X_vnir(valid,:);
X_swir = X_swir(valid,:);
Y = Y(valid);
nFeat = size(X_raw, 2);


%Converting to table
predNames = compose("WL_%.1f", WL_combined(:));
predNames_vnir = compose("VNIR_%.1f", WL_vnir(:));
predNames_swir = compose("SWIR_%.1f", WL_swir(9:end));

Tbl_raw = array2table(X_raw, "VariableNames", cellstr(predNames));
Tbl_raw.Label = Y;

Tbl_vnir = array2table(X_vnir, "VariableNames", cellstr(predNames_vnir));
Tbl_vnir.Label = Y;

Tbl_swir = array2table(X_swir, "VariableNames", cellstr(predNames_swir));
Tbl_swir.Label = Y;

%% Preprocessing variants !! for the cases where I want to test the set with the independent image;
% Preprocessing variants

% Training set only
X_snv_vnir = SNV(X_vnir);
X_msc_vnir = MSC_reworked(X_vnir, mean(X_vnir,1));

X_snv_swir = SNV(X_swir);
X_msc_swir = MSC_reworked(X_swir, mean(X_swir,1));

X_snv_combined = SNV(X_raw);
X_msc_combined = MSC_reworked(X_raw, mean(X_raw,1));




% Testing set only (independent image)
% Raw
X_vnir_mix = Tbl_vnir_mean_mix{:,1:end-1};
X_swir_mix = Tbl_swir_mean_mix{:,9:end-1};
X_combined_mix = Tbl_combined_mix{:,1:end-1};

% Pretreatment for VNIR
X_snv_mix_VNIR = SNV(X_vnir_mix);
X_msc_mix_VNIR = MSC_reworked(X_vnir_mix, mean(X_vnir,1));

% Pretreatment for SWIR
X_snv_mix_SWIR = SNV(X_swir_mix);
X_msc_mix_SWIR = MSC_reworked(X_swir_mix, mean(X_swir,1));

% Pretreatment for combined
X_snv_mix_combined = SNV(X_combined_mix);
X_msc_mix_combined = MSC_reworked(X_combined_mix, mean(X_raw,1));


%Fixing the issue with the labels
Y_mix = Tbl_combined_mix.Label;
if ~iscategorical(Y_mix)
    Y_mix = categorical(Y_mix);
end
Y_mix = renamecats(Y_mix, "Healthy", "HLTHY");
Y_mix = categorical(string(Y_mix), categories(Y));
missingMixLabels = isundefined(Y_mix);
assert(~any(missingMixLabels), ...
    "Independent test data contains classes not present in training data: %s", ...
    strjoin(unique(string(Tbl_combined_mix.Label(missingMixLabels))), ", "));

predNames_vnir_mix = predNames_vnir;
predNames_swir_mix = predNames_swir;
predNames_combined_mix = predNames;


%% Creating tables for Classification Learner
% Training set only

%combined dataset
Tbl_raw = array2table(X_raw, "VariableNames", cellstr(predNames));
Tbl_raw.Label = Y;

Tbl_vnir_raw = array2table(X_vnir, "VariableNames", cellstr(predNames_vnir));
Tbl_vnir_raw.Label = Y;

Tbl_swir_raw = array2table(X_swir, "VariableNames", cellstr(predNames_swir));
Tbl_swir_raw.Label = Y;

Tbl_snv_vnir = array2table(X_snv_vnir, "VariableNames", cellstr(predNames_vnir)); 
Tbl_snv_vnir.Label = Y;

Tbl_msc_vnir = array2table(X_msc_vnir, "VariableNames", cellstr(predNames_vnir)); 
Tbl_msc_vnir.Label = Y;

Tbl_snv_swir = array2table(X_snv_swir, "VariableNames", cellstr(predNames_swir)); 
Tbl_snv_swir.Label = Y;

Tbl_msc_swir = array2table(X_msc_swir, "VariableNames", cellstr(predNames_swir)); 
Tbl_msc_swir.Label = Y;

Tbl_snv_combined = array2table(X_snv_combined, "VariableNames", cellstr(predNames));
Tbl_snv_combined.Label = Y;

Tbl_msc_combined = array2table(X_msc_combined, "VariableNames", cellstr(predNames));
Tbl_msc_combined.Label = Y;




% Testing set only (independent image)
Tbl_raw_mix = array2table(X_combined_mix, "VariableNames", cellstr(predNames_combined_mix));
Tbl_raw_mix.Label = Y_mix;

Tbl_vnir_mix_raw = array2table(X_vnir_mix, "VariableNames", cellstr(predNames_vnir_mix));
Tbl_vnir_mix_raw.Label = Y_mix;

Tbl_swir_mix_raw = array2table(X_swir_mix, "VariableNames", cellstr(predNames_swir_mix));
Tbl_swir_mix_raw.Label = Y_mix;


Tbl_snv_mix_VNIR = array2table(X_snv_mix_VNIR, "VariableNames", cellstr(predNames_vnir_mix));
Tbl_snv_mix_VNIR.Label = Y_mix;

Tbl_msc_mix_VNIR = array2table(X_msc_mix_VNIR, "VariableNames", cellstr(predNames_vnir_mix));
Tbl_msc_mix_VNIR.Label = Y_mix;

Tbl_snv_mix_SWIR = array2table(X_snv_mix_SWIR, "VariableNames", cellstr(predNames_swir_mix));
Tbl_snv_mix_SWIR.Label = Y_mix;

Tbl_msc_mix_SWIR = array2table(X_msc_mix_SWIR, "VariableNames", cellstr(predNames_swir_mix));
Tbl_msc_mix_SWIR.Label = Y_mix;

Tbl_snv_mix_combined = array2table(X_snv_mix_combined, "VariableNames", cellstr(predNames_combined_mix));
Tbl_snv_mix_combined.Label = Y_mix;

Tbl_msc_mix_combined = array2table(X_msc_mix_combined, "VariableNames", cellstr(predNames_combined_mix));
Tbl_msc_mix_combined.Label = Y_mix;







% plotting the mean spectra for the msc corrected dataset 
Tbl_msc_combined = array2table(X_msc_combined, "VariableNames", cellstr(predNames));
Tbl_msc_combined.Label = Y;

plotClassMeanSpectra(Tbl_msc_combined, predNames, "MSC Corrected Spectra (Training Set)");



figure, 
plotClassMeanSpectra(Tbl_msc_mix_combined, predNames_combined_mix, "MSC Corrected Spectra (Independent Test Set)");
























% %% Creating stratified training and test sets
% 
% testFraction = 0.30;
% rng(1); % Reproducible split
% 
% cv = cvpartition(Y, "HoldOut", testFraction);
% idxTrain = training(cv);
% idxTest = test(cv);
% 
% Tbl_raw_train = Tbl_raw(idxTrain,:);
% Tbl_raw_test  = Tbl_raw(idxTest,:);
% 
% disp("Training class counts:");
% disp(groupcounts(Tbl_raw_train.Label));
% 
% % Pretreatment on training set only
% X_snv_train = SNV(Tbl_raw_train{:,1:nFeat});
% X_msc_train = MSC_reworked(Tbl_raw_train{:,1:nFeat}, mean(Tbl_raw_train{:,1:nFeat},1));
% 
% Tbl_snv_train = array2table(X_snv_train, "VariableNames", cellstr(predNames));
% Tbl_snv_train.Label = Tbl_raw_train.Label;
% 
% Tbl_msc_train = array2table(X_msc_train, "VariableNames", cellstr(predNames));
% Tbl_msc_train.Label = Tbl_raw_train.Label;
% 
% 
% 
% %% Applying the same pretreatment to test set (using training set parameters)
% X_snv_test = SNV(Tbl_raw_test{:,1:nFeat});
% X_msc_test = MSC_reworked(Tbl_raw_test{:,1:nFeat}, mean(Tbl_raw_train{:,1:nFeat},1));
% 
% 
% Tbl_snv_test = array2table(X_snv_test, "VariableNames", cellstr(predNames));
% Tbl_snv_test.Label = Tbl_raw_test.Label;
% 
% 
% Tbl_msc_test = array2table(X_msc_test, "VariableNames", cellstr(predNames));
% Tbl_msc_test.Label = Tbl_raw_test.Label;
% 
% 
% 
% %% Plot class-wise mean spectra for RAW, SNV, MSC
% 
% classes = categories(Y);
% nClasses = numel(classes);
% 
% % X-axis: feature index (replace with wavelength vector if you have one)
% x = WL_combined(:);
% 
% % Preallocate class means
% mean_raw = nan(nClasses, nFeat);
% mean_snv = nan(nClasses, nFeat);
% mean_msc = nan(nClasses, nFeat);
% 
% for k = 1:nClasses
%     idx = (Y == classes{k});
%     mean_raw(k,:) = mean(X_raw(idx,:), 1, 'omitnan');
%     mean_snv(k,:) = mean(X_snv_combined(idx,:), 1, 'omitnan');
%     mean_msc(k,:) = mean(X_msc_combined(idx,:), 1, 'omitnan');
% end
% 
% figure('Color','w');
% tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');
% 
% % RAW
% nexttile; hold on;
% for k = 1:nClasses
%     plot(x, mean_raw(k,:), 'LineWidth', 1.6);
% end
% hold off; grid on;
% title('RAW');
% xlabel('Wavelength (nm)');
% ylabel('Reflectance');
% legend(classes, 'Location', 'best');
% 
% % SNV
% nexttile; hold on;
% for k = 1:nClasses
%     plot(x, mean_snv(k,:), 'LineWidth', 1.6);
% end
% hold off; grid on;
% title('SNV');
% xlabel('Wavelength (nm)');
% ylabel('SNV Value');
% legend(classes, 'Location', 'best');
% 
% % MSC
% nexttile; hold on;
% for k = 1:nClasses
%     plot(x, mean_msc(k,:), 'LineWidth', 1.6);
% end
% hold off; grid on;
% title('MSC');
% xlabel('Wavelength (nm)');
% ylabel('MSC Corrected Value');
% legend(classes, 'Location', 'best');
% 
% sgtitle('Class-wise Mean Spectra Across Pretreatments');
% 
% 
% %% Plot class-wise mean spectra (+/- SD) for RAW, SNV, MSC
% 
% classes = categories(Y);
% nClasses = numel(classes);
% x = WL_combined(:);   % wavelength axis
% 
% % Preallocate class means and SDs
% mean_raw = nan(nClasses, nFeat);
% mean_snv = nan(nClasses, nFeat);
% mean_msc = nan(nClasses, nFeat);
% 
% std_raw  = nan(nClasses, nFeat);
% std_snv  = nan(nClasses, nFeat);
% std_msc  = nan(nClasses, nFeat);
% 
% for k = 1:nClasses
%     idx = (Y == classes{k});
% 
%     mean_raw(k,:) = mean(X_raw(idx,:), 1, 'omitnan');
%     mean_snv(k,:) = mean(X_snv_combined(idx,:), 1, 'omitnan');
%     mean_msc(k,:) = mean(X_msc_combined(idx,:), 1, 'omitnan');
% 
%     std_raw(k,:)  = std(X_raw(idx,:), 0, 1, 'omitnan');
%     std_snv(k,:)  = std(X_snv_combined(idx,:), 0, 1, 'omitnan');
%     std_msc(k,:)  = std(X_msc_combined(idx,:), 0, 1, 'omitnan');
% end
% 
% figure('Color','w');
% tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');
% 
% clr = lines(nClasses);
% alphaBand = 0.15;
% 
% % RAW
% nexttile; hold on;
% for k = 1:nClasses
%     mu = mean_raw(k,:);
%     sd = std_raw(k,:);
%     c  = clr(k,:);
% 
%     fill([x; flipud(x)]', [mu-sd, fliplr(mu+sd)], c, ...
%         'FaceAlpha', alphaBand, 'EdgeColor', 'none');
%     plot(x, mu, 'LineWidth', 1.8, 'Color', c);
% end
% hold off; grid on;
% title('RAW');
% xlabel('Wavelength (nm)');
% ylabel('Reflectance');
% legend(classes, 'Location', 'best');
% 
% % SNV
% nexttile; hold on;
% for k = 1:nClasses
%     mu = mean_snv(k,:);
%     sd = std_snv(k,:);
%     c  = clr(k,:);
% 
%     fill([x; flipud(x)]', [mu-sd, fliplr(mu+sd)], c, ...
%         'FaceAlpha', alphaBand, 'EdgeColor', 'none');
%     plot(x, mu, 'LineWidth', 1.8, 'Color', c);
% end
% hold off; grid on;
% title('SNV');
% xlabel('Wavelength (nm)');
% ylabel('SNV Value');
% legend(classes, 'Location', 'best');
% 
% % MSC
% nexttile; hold on;
% for k = 1:nClasses
%     mu = mean_msc(k,:);
%     sd = std_msc(k,:);
%     c  = clr(k,:);
% 
%     fill([x; flipud(x)]', [mu-sd, fliplr(mu+sd)], c, ...
%         'FaceAlpha', alphaBand, 'EdgeColor', 'none');
%     plot(x, mu, 'LineWidth', 1.8, 'Color', c);
% end
% hold off; grid on;
% title('MSC');
% xlabel('Wavelength (nm)');
% ylabel('MSC Corrected Value');
% legend(classes, 'Location', 'best');
% 
% sgtitle('Class-wise Mean Spectra Across Pretreatments (Mean +/- 1 SD)');




