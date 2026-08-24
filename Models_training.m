% cd('')

clear variables
close all 

%% Load the VNIR and SWIR summary tables
% sourcedir = "";

srcFolder = uigetdir(sourcedir,'source folder');
files = dir(fullfile(srcFolder, '*.mat'));

for i = 1:numel(files)
    filePath = fullfile(srcFolder, files(i).name);
    load(filePath);
end 



%% choosing the complexity
%reduced complexity for horizontal concatenation + the specific orientation 
% isTarget = @(name) (contains(name, '_BEM_') | contains(name, '_HLTHY_')); ...
% % & contains(name, '_up_');

% trainMask = isTarget(Tbl_combined_samples.Sample_Name_A) | isTarget(Tbl_combined_samples.Sample_Name_B);

% testMask = isTarget(Tbl_combined_mix.Sample_Name_A) | isTarget(Tbl_combined_mix.Sample_Name_B);

% table_train = Tbl_combined_samples(trainMask, :);
% table_test  = Tbl_combined_mix(testMask, :);

% skipping the step if we go with the original complexity




%% in case the horizontal concatenation is used
% up-side orientation
% table_train = Tbl_combined_samples(contains(Tbl_combined_samples.Sample_Name_A, '_up_')|contains(Tbl_combined_samples.Sample_Name_B, '_up_'), :);
% table_test = Tbl_combined_mix(contains(Tbl_combined_mix.Sample_Name_A, '_up_')|contains(Tbl_combined_mix.Sample_Name_B, '_up_'), :);



% down-side orientation
% table_train = Tbl_combined_samples(contains(Tbl_combined_samples.Sample_Name_A, '_down_')|contains(Tbl_combined_samples.Sample_Name_B, '_down_'), :);
% table_test = Tbl_combined_mix(contains(Tbl_combined_mix.Sample_Name_A, '_down_')|contains(Tbl_combined_mix.Sample_Name_B, '_down_'), :);



% vertical concatenation
table_train = Tbl_combined_samples;
table_test = Tbl_combined_mix;
%% Side orientation 

%up-side orientation 
% table_train = SWIR_Summary(contains(SWIR_Summary.Sample_Name,'up'),:);
% table_test = MIX_SWIR_Summary(contains(MIX_SWIR_Summary.Sample_Name,'up'),:);


% down-side orientation 
% table_train = SWIR_Summary(contains(SWIR_Summary.Sample_Name,'down'),:);
% table_test = MIX_SWIR_Summary(contains(MIX_SWIR_Summary.Sample_Name,'down'),:);


% vertical concatenation
% table_train = SWIR_Summary;
% table_test = MIX_SWIR_Summary;


%% Extracting the labels and preparing the tables
% IN CASE HORIZONTAL CONCATENATION IS USED
train_parts = split(string(table_train.Sample_Name_A), "_");
test_parts  = split(string(table_test.Sample_Name_A), "_");


% train_parts = split(string(table_train.Sample_Name), "_");
% test_parts  = split(string(table_test.Sample_Name), "_");

%% Separating labels
train_labels = train_parts(:,2);
test_labels  = test_parts(:,2);

% train_labels(ismember(train_labels, ["PEN"])) = "PEM";
% test_labels(ismember(test_labels, ["PEN"])) = "PEM";


% BAD LABELLING
train_labels(ismember(train_labels, ["BEM", "BEN", "PEN"])) = "BAD";
test_labels(ismember(test_labels, ["BEM", "BEN", "PEN"])) = "BAD";


class_names = unique(train_labels, "stable");

[is_train, Y_train] = ismember(train_labels, class_names);
[is_test,  Y_test]  = ismember(test_labels, class_names);


%%  Training
% IN CASE HORIZONTAL CONCATENATION IS USED
X_train = cell2mat(table_train.Combined_Mean_Spectrum);


% X_train = cell2mat(table_train.Mean_Spectrum);
%%  Test
% IN CASE HORIZONTAL CONCATENATION IS USED
X_test = cell2mat(table_test.Combined_Mean_Spectrum);


% X_test = cell2mat(table_test.Mean_Spectrum);


%% Plotting the spectra for the training set
figure;
bandIndex = table_train.Combined_Wavelengths{1};

figure('Color', 'w', 'Name', 'Full VNIR-SWIR Spectra for vertical orientation');
hold on;
sampleLines = plot(bandIndex, X_train', 'LineWidth', 1, 'HandleVisibility', 'off');
meanLine = plot(bandIndex, mean(X_train, 1), 'k--', 'LineWidth', 2, ...
    'HandleVisibility', 'off');
hold off;
xlabel('Wavelength (nm)');
ylabel('Reflectance (a.u.)');
title('VNIR-SWIR Spectra for Samples (vertical orientation)');
grid on;
legend([sampleLines(1), meanLine], {'Samples', 'Mean'}, 'Location', 'best');

%% RUN the PLSDA

% PLS-DA: train on samples, evaluate on independent mix test set
[output_PLSDA, ~, ~] = PLSDA_new(X_train, Y_train, X_test, Y_test, 20, 100);

%% Calibration results
output_PLSDA{1,1}.CC_cal(6)
output_PLSDA{1,4}.CC_cal(6)
output_PLSDA{1,6}.CC_cal(8)
%% Test results
% RAW, 2D and SNV+2D

output_PLSDA{1,1}.SS{1,6}.CC
output_PLSDA{1,1}.SS{1,6}.F1

output_PLSDA{1,4}.SS{1,6}.CC
output_PLSDA{1,4}.SS{1,6}.F1

output_PLSDA{1,6}.SS{1,8}.CC
output_PLSDA{1,6}.SS{1,8}.F1

%% displaying the confusion matrix with labels 
% After PLSDA runs
Pred_Ytest = output_PLSDA{1,6}.Pred_class_t(:,6);
class_names_plot = unique(train_labels, "stable");

C = confusionmat(Y_test, Pred_Ytest);
figure;
confusionchart(C, class_names_plot);
title('Test Set Confusion Matrix');


%% doing the pretreatment for classification learner app 
% After the PLSDA call, apply the same pretreatment manually:
X_train_2D = savgol(X_train, 15, 3, 2);
X_test_2D = savgol(X_test,  15, 3, 2);


X_train_SNV_SG2D = savgol(SNV(X_train), 15, 3, 2);
X_test_SNV_SG2D = savgol(SNV(X_test),  15, 3, 2);


%% Mean overall spectra
% 
% 2D pretreated spectra for the training set
figure;
figure('Color', 'w', 'Name', '2D pretreated Spectra for up-side orientation');
plot(bandIndex, X_test_2D', 'LineWidth', 1);
hold on;
hold off;
xlabel('Wavelength (nm)');
ylabel('Second derivative reflectance (a.u.)');
title('2D pretreated VNIR-SWIR Spectra for Samples (up-side orientation)');
grid on;


% SNV+2D pretreated spectra for the training set
figure;
figure('Color', 'w', 'Name', 'SNV+2D pretreated Spectra for up-side orientation');
plot(bandIndex, X_test_SNV_SG2D', 'LineWidth', 1);
hold on;
hold off;
xlabel('Wavelength (nm)');
ylabel('Second derivative of SNV-transformed reflectance (a.u.)');
title('SNV+2D pretreated VNIR-SWIR Spectra for Samples (up-side orientation)');
grid on;


%% Class-wise mean spectra 
% RAW
figure('Color', 'w', 'Name', 'Mean Spectra per Class');
hold on;
for c = 1:numel(class_names)
    classMask = train_labels == class_names(c);
    meanSpec = mean(X_train(classMask, :), 1);
    plot(bandIndex, meanSpec, 'LineWidth', 2, 'DisplayName', class_names(c));
end
hold off;
legend;
xlabel('Wavelength (nm)');
ylabel('Reflectance (a.u.)');
title('Mean Spectra per Class (up-side orientation)');
grid on;

% 2D 
figure('Color', 'w', 'Name', 'Mean Spectra per Class');
hold on;
for c = 1:numel(class_names)
    classMask = train_labels == class_names(c);
    meanSpec = mean(X_train_2D(classMask, :), 1);
    plot(bandIndex, meanSpec, 'LineWidth', 2, 'DisplayName', class_names(c));
end
hold off;
legend;
xlabel('Wavelength (nm)');
ylabel('Second derivative of reflectance (a.u.)');
title(' 2D pretreated VNIR-SWIR Spectra for Samples per Class (up-side orientation)');
grid on;

% 2D SNV 
figure('Color', 'w', 'Name', 'Mean Spectra per Class');
hold on;
for c = 1:numel(class_names)
    classMask = train_labels == class_names(c);
    meanSpec = mean(X_train_SNV_SG2D(classMask, :), 1);
    plot(bandIndex, meanSpec, 'LineWidth', 2, 'DisplayName', class_names(c));
end
hold off;
legend;
xlabel('Wavelength (nm)');
ylabel('Second derivative of SNV-transformed reflectance (a.u.)');
title('SNV+2D pretreated VNIR-SWIR Spectra for Samples per Class (up-side orientation)');
grid on;

















