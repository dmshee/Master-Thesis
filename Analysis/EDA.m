cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Analysis');
load('MIX_classification_inputs.mat');

%% Analsying the aborbance spectra for the MIX samples

% Extract wavelength axis and spectra
wavelengths = MIX_sampleMean(1).wavelengths(:);
spectra = vertcat(MIX_sampleMean.meanSpectrum);

% Extract category labels
classNames = string({MIX_sampleMean.className});

% Normalize absorbance data per sample spectrum

spectraNorm = (spectra - mean(spectra, 2, "omitnan")) ./ std(spectra, 0, 2, "omitnan");

% Plot one combined mean spectrum per category
categories = unique(classNames, "stable");

for i = 1:numel(categories)
    category = categories(i);
    idx = classNames == category;

    categoryMeanSpectrum = mean(spectraNorm(idx, :), 1, "omitnan");

    figure;
    plot(wavelengths, categoryMeanSpectrum, "LineWidth", 2);
    grid on;

    xlabel("Wavelength");
    ylabel("Normalized absorbance");
    title("Combined mean spectrum - " + category);
end

figure;
plot(wavelengths, spectraNorm(idx, :)', "Color", [0.75 0.75 0.75]);
hold on;
plot(wavelengths, categoryMeanSpectrum, "k", "LineWidth", 2);
grid on;

xlabel("Wavelength");
ylabel("Normalized absorbance");
title("Combined mean spectrum - " + category);
legend("Samples", "Category mean");




%% Analysing the data with PCA
X = Tbl_combined_mix{:, startsWith(Tbl_combined_mix.Properties.VariableNames, 'WL_')};

groups =Tbl_combined_mix.Label;

[L,S, Ev] = PCA(X);

explained_variance = Ev / sum(Ev) * 100;

%PC plots 

pcPairs = [1 2; 1 3; 2 3];

for i = 1:size(pcPairs, 1)

    pcX = pcPairs(i, 1);
    pcY = pcPairs(i, 2);

    figure;
    gscatter(S(:,pcX), S(:,pcY), groups, [], '.', 22);

    xlabel(sprintf('PC%d (%.2f%%)', pcX, explained_variance(pcX)));
    ylabel(sprintf('PC%d (%.2f%%)', pcY, explained_variance(pcY)));
    title(sprintf('PCA Scatter Plot: PC%d vs PC%d', pcX, pcY));

    grid on;
    legend('Location', 'best');

end


%% Seeing the scatter for the pretreated data

SNV_X = SNV(X);

[L_SNV,S_SNV, Ev_SNV] = PCA(SNV_X);

for i = 1:size(pcPairs, 1)

    pcX = pcPairs(i, 1);
    pcY = pcPairs(i, 2);

    figure;
    gscatter(S_SNV(:,pcX), S_SNV(:,pcY), groups, [], '.', 22);

    xlabel(sprintf('PC%d (%.2f%%)', pcX, explained_variance(pcX)));
    ylabel(sprintf('PC%d (%.2f%%)', pcY, explained_variance(pcY)));
    title(sprintf('PCA Scatter Plot SNV: PC%d vs PC%d', pcX, pcY));

    grid on;
    legend('Location', 'best');

end


%% Extracting the mean spectra for the training class 
% Choosing the spectral range for analysis: VNIR, SWIR, or combined
X_vnir = cell2mat(T_merged.VNIR_Mean_Spectrum);
X_swir = cell2mat(T_merged.SWIR_Mean_Spectrum);
X_combined = cell2mat(T_merged.Combined_Spectrum);

class_text = regexp(T_merged.Sample_ID, '^[^_]+', 'match', 'once');
class_cat = categorical(class_text);
class_names = categories(class_cat);
Y = double(class_cat);