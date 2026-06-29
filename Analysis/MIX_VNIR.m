clear variables
close all


cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Analysis')

%% Displaying the mean image across all wavelength bands for the Mix reflectance data
load('MIX_VNIR.mat')

mean_image = mean(im_reflectance, 3);
figure, imshow(mean_image, [])
title('Mean Image Across All Wavelength Bands for Mix Reflectance VNIR')

figure, imshow(im_reflectance(:,:,172), []);
title('Image at 172nd Wavelength Band (Mix Reflectance VNIR)')


%% Alternative approach to create a pseudo-RGB image by manually selecting three bands
% Choose three bands manually
R_band = 134;
G_band = 145;
B_band = 151;

% Extract the three bands
R = im_reflectance(:, :, R_band);
G = im_reflectance(:, :, G_band);
B = im_reflectance(:, :, B_band);

% Combine into pseudo-RGB image
pseudoRGB = cat(3, R, G, B);

% Rescale values to [0, 1] for display
pseudoRGB = rescale(pseudoRGB);

% Display
figure;
imshow(pseudoRGB);
title('Pseudo-RGB Image from Selected Hyperspectral Bands for Mix Reflectance VNIR');





%% Creating mask through the segmentation of the samples
% taking a look at the spectra for the image

% standard deviation image

Xuf = unfold(im_reflectance);     
im_std = std(Xuf,0,2);            
im_sd = reshape(im_std, size(im_reflectance,1), size(im_reflectance,2));

%Finding the threshold and making the mask 
figure
subplot (1,2,1)
imshow(im_sd,[]);

subplot(1,2,2)
histogram(im_sd,100);
title('Pixel Intensities of standard deviation image')

% Applying the mask to the mean image
mask1 = im_sd>1.56;
figure, imshow(mask1.*mean_image,[])
title('Mean Image with Mask Applied')

% Displaying the mask and filling in the holes
figure,
imshow(mask1,[]);

mask2 = imfill(mask1,"holes");
figure,
imshow(mask2,[]);

%Labelling the samples
L= bwlabel(mask2);
figure, 
imshow(L,[]);

colormap('jet')


%% Commiting to alternative approach 
x = zeros(1,25); 
y = zeros(1,25);  

for i = 1:25
    fprintf('Select the central region of sample %d.\n',i);
    [ xi, yi, ~] = impixel(im_reflectance(:,:,172),[]);
    %storing columns (x coordinates) 
    x(i) = round(xi(1)); 
    % Storing rows (y coordinates) 
    y(i) = round(yi(1));  
end



% Creating a matrix to store sample labels
L_new = zeros(size(L));  

for i = 1:25
    % Getting the old label at the selected coordinates
    L_old = L(y(i), x(i));
    % Replacing all pixels with the old label with the new label i 
    L_new(L == L_old) = i; 
end
unique_labels = unique(L_new);
disp('Unique labels in L_new:');
disp(unique_labels);


%% Creating the class names 
classNames = {'Healthy', 'BEM','BEN','PEN','ODEFF'};
classSamples = [1,7,13,19,25; 2,3,8,14,20; 6,11,12,18,24; 16,17,21,22,23; 4,5,9,10,15];


%% Building the class map from the sample labels

% Creating a vector to store class labels 
sample_to_class = zeros(1,25);  
for classID = 1:length(classNames)
    % Assigning class IDs 
    sample_to_class(classSamples(classID, :)) = classID;  
end

%Building class map image from sample labels
class_map = zeros(size(L_new));
for s = 1:25
    class_map(L_new == s) = sample_to_class(s);
end

% Displaying the class map
figure,
imshow(class_map, []);
colormap('jet');
colorbar('Ticks', 1:length(classNames), 'TickLabels', classNames);
title('Class Map of Samples')


%% Creating the structure with the mean spectra for each sample

sampleMap = L_new; 
classMap = class_map; 


sampleMean_VNIR = struct();

for sampleID = 1:25
    classID = sample_to_class(sampleID);
    mask = sampleMap == sampleID;
    mask_lin = mask(:);
    sampleSpectra = Xuf(mask_lin, :);
    

    % Calculating the mean spectrum for the sample (mean across all pixels in the sample)
    meanSpectrum = zeros(1, numel(WL));

    for b = 1:numel(WL)
        bandImage = im_reflectance(:, :, b);
        meanSpectrum(b) = mean(bandImage(mask), "omitnan");
    end


    sampleMean_VNIR(sampleID).sampleID = sampleID;
    sampleMean_VNIR(sampleID).sampleName = sprintf('%s_%02d', classNames{classID}, sampleID);
    sampleMean_VNIR(sampleID).classID = classID;
    sampleMean_VNIR(sampleID).className = classNames{classID};
    sampleMean_VNIR(sampleID).wavelengths = WL;
    sampleMean_VNIR(sampleID).meanSpectrum = meanSpectrum;
    sampleMean_VNIR(sampleID).unfoldedSpectra = sampleSpectra;
    sampleMean_VNIR(sampleID).nPixels = size(sampleSpectra, 1);
end


% saving the sample mean spectra and class map for later use
save('MIX_VNIR_sampleMean.mat', 'sampleMean_VNIR', 'classMap', 'WL');













































%% Alternative stash 

% %% Preprocessing
% % Apply SNV
% sample_mean_spectra_snv = SNV(sample_mean_spectra);

% % Apply MSC (needs a target mean - use the global mean of all spectra)
% target_mean = mean(sample_mean_spectra, 1);
% sample_mean_spectra_msc = MSC_reworked(sample_mean_spectra, target_mean);

% % Plot SNV-preprocessed class means
% figure
% hold on
% for classID = 1:numel(classNames)
%     sampleIDs = classSamples(classID, :);
%     class_mean_snv = mean(sample_mean_spectra_snv(sampleIDs, :), 1);
%     plot(WL, class_mean_snv, 'LineWidth', 2)
% end
% hold off
% grid on
% xlabel('Wavelength (nm)')
% ylabel('SNV-Normalized Reflectance')
% title('Mean Spectral Profiles by Class (SNV)')
% legend(classNames, 'Location', 'best')

% % Plot MSC-preprocessed class means
% figure
% hold on
% for classID = 1:numel(classNames)
%     sampleIDs = classSamples(classID, :);
%     class_mean_msc = mean(sample_mean_spectra_msc(sampleIDs, :), 1);
%     plot(WL, class_mean_msc, 'LineWidth', 2)
% end
% hold off
% grid on
% xlabel('Wavelength (nm)')
% ylabel('MSC-Normalized Reflectance')
% title('Mean Spectral Profiles by Class (MSC)')
% legend(classNames, 'Location', 'best')


% %% Doing PCA on the masked image 
% % PCA comparison: untreated vs SNV vs MSC
% X_raw = sample_mean_spectra;      
% X_snv = SNV(X_raw);
% X_msc = MSC_reworked(X_raw, mean(X_raw, 1));

% [L_raw, S_raw, Ev_raw] = PCA(X_raw);
% [L_snv, S_snv, Ev_snv] = PCA(X_snv);
% [L_msc, S_msc, Ev_msc] = PCA(X_msc);

% % Variance explained (untreated, SNV, and MSC)
% variance_raw = 100 * Ev_raw / sum(Ev_raw);
% variance_snv = 100 * Ev_snv / sum(Ev_snv);
% variance_msc = 100 * Ev_msc / sum(Ev_msc);

% nPC = min(5, numel(variance_raw));   % safe if fewer than 5 PCs

% figure
% plot(1:nPC, variance_raw(1:nPC), '--o', 'LineWidth', 1.5)
% hold on
% plot(1:nPC, variance_snv(1:nPC), '--s', 'LineWidth', 1.5)
% plot(1:nPC, variance_msc(1:nPC), '--^', 'LineWidth', 1.5)
% hold off
% grid on
% xlabel('PC')
% ylabel('% variance explained')
% title('Variance Explained: Untreated vs SNV vs MSC')
% legend('Untreated', 'SNV', 'MSC', 'Location', 'best')


% %% Plotting the PCA applied on the masked images
% [nRows, nCols, nBands] = size(im_reflectance);

% % Unfold full image to pixels x bands
% X_all = unfold(im_reflectance);

% % Use labeled sample area as mask (same region you used for class mapping)
% mask_pixels = L_new(:) > 0;
% X_masked = X_all(mask_pixels, :);

% % PCA on untreated masked spectra
% [L_raw_px, S_raw_px, Ev_raw_px] = PCA(X_masked);

% % PCA on SNV-preprocessed masked spectra
% X_masked_snv = SNV(X_masked);
% [L_snv_px, S_snv_px, Ev_snv_px] = PCA(X_masked_snv);

% % PCA on MSC-preprocessed masked spectra
% X_masked_msc = MSC_reworked(X_masked, mean(X_masked, 1));
% [L_msc_px, S_msc_px, Ev_msc_px] = PCA(X_masked_msc);

% % Rebuild PC1-3 score images to full image size (NaN background)
% S_raw_full = nan(nRows*nCols, 3);
% S_raw_full(mask_pixels, :) = S_raw_px(:, 1:3);
% S_raw_im = reshape(S_raw_full, nRows, nCols, 3);

% S_snv_full = nan(nRows*nCols, 3);
% S_snv_full(mask_pixels, :) = S_snv_px(:, 1:3);
% S_snv_im = reshape(S_snv_full, nRows, nCols, 3);

% S_msc_full = nan(nRows*nCols, 3);
% S_msc_full(mask_pixels, :) = S_msc_px(:, 1:3);
% S_msc_im = reshape(S_msc_full, nRows, nCols, 3);

% % Plot untreated, SNV, and MSC, first 3 PCs
% figure('Position', [100 100 1400 900])
% colormap(jet)   
% for pc = 1:3
%     subplot(3,3,pc)
%     imagesc(S_raw_im(:,:,pc))
%     axis image off
%     colorbar
%     title(sprintf('Untreated PC%d', pc))

%     subplot(3,3,pc+3)
%     imagesc(S_snv_im(:,:,pc))
%     axis image off
%     colorbar
%     title(sprintf('SNV PC%d', pc))
    
%     subplot(3,3,pc+6)
%     imagesc(S_msc_im(:,:,pc))
%     axis image off
%     colorbar
%     title(sprintf('MSC PC%d', pc))
% end
% sgtitle('First 3 PCA Score Images on Masked Region (Untreated vs SNV vs MSC)')

% %% Running K-means on the masked spectra (untreated vs SNV vs MSC)

% % Number of clusters (use known number of classes)
% k = numel(classNames);

% % Reproducibility
% rng(1);

% % K-means on untreated masked spectra
% [idx_raw, C_raw] = kmeans(X_masked, k, ...
%     'Distance', 'sqeuclidean', ...
%     'Replicates', 20, ...
%     'MaxIter', 500);

% % K-means on SNV-preprocessed masked spectra
% [idx_snv, C_snv] = kmeans(X_masked_snv, k, ...
%     'Distance', 'sqeuclidean', ...
%     'Replicates', 20, ...
%     'MaxIter', 500);

% % K-means on MSC-preprocessed masked spectra
% [idx_msc, C_msc] = kmeans(X_masked_msc, k, ...
%     'Distance', 'sqeuclidean', ...
%     'Replicates', 20, ...
%     'MaxIter', 500);

% % Map cluster labels back to full image (0 = background)
% km_raw_full = zeros(nRows*nCols, 1);
% km_raw_full(mask_pixels) = idx_raw;
% km_raw_im = reshape(km_raw_full, nRows, nCols);

% km_snv_full = zeros(nRows*nCols, 1);
% km_snv_full(mask_pixels) = idx_snv;
% km_snv_im = reshape(km_snv_full, nRows, nCols);

% km_msc_full = zeros(nRows*nCols, 1);
% km_msc_full(mask_pixels) = idx_msc;
% km_msc_im = reshape(km_msc_full, nRows, nCols);

% % Display cluster maps
% figure('Position', [100 100 1500 450])

% subplot(1,3,1)
% imagesc(km_raw_im)
% axis image off
% title('K-means Clusters (Untreated)')
% colormap(parula(k))
% cb1 = colorbar;
% cb1.Ticks = 1:k;
% cb1.TickLabels = compose('Cluster %d', 1:k);

% subplot(1,3,2)
% imagesc(km_snv_im)
% axis image off
% title('K-means Clusters (SNV)')
% colormap(parula(k))
% cb2 = colorbar;
% cb2.Ticks = 1:k;
% cb2.TickLabels = compose('Cluster %d', 1:k);

% subplot(1,3,3)
% imagesc(km_msc_im)
% axis image off
% title('K-means Clusters (MSC)')
% colormap(parula(k))
% cb3 = colorbar;
% cb3.Ticks = 1:k;
% cb3.TickLabels = compose('Cluster %d', 1:k);

% % Plot mean spectrum per cluster (untreated)
% figure
% hold on
% for c = 1:k
%     plot(WL, mean(X_masked(idx_raw == c, :), 1), 'LineWidth', 2)
% end
% hold off
% grid on
% xlabel('Wavelength (nm)')
% ylabel('Mean Reflectance')
% title('Cluster Mean Spectra (Untreated)')
% legend(compose('Cluster %d', 1:k), 'Location', 'best')

% % Plot mean spectrum per cluster (SNV)
% figure
% hold on
% for c = 1:k
%     plot(WL, mean(X_masked_snv(idx_snv == c, :), 1), 'LineWidth', 2)
% end
% hold off
% grid on
% xlabel('Wavelength (nm)')
% ylabel('SNV Reflectance')
% title('Cluster Mean Spectra (SNV)')
% legend(compose('Cluster %d', 1:k), 'Location', 'best')

% % Plot mean spectrum per cluster (MSC)
% figure
% hold on
% for c = 1:k
%     plot(WL, mean(X_masked_msc(idx_msc == c, :), 1), 'LineWidth', 2)
% end

% % hold off
% % grid on
% % xlabel('Wavelength (nm)')
% % ylabel('MSC Reflectance')
% % title('Cluster Mean Spectra (MSC)')
% % legend(compose('Cluster %d', 1:k), 'Location', 'best')
