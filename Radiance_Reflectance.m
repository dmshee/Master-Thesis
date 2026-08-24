% cd('');

srcFolder = uigetdir('','source folder');
dstFolder = uigetdir('', 'destination folder');

files = dir(fullfile(srcFolder, '*.img'));

for k = 1:numel(files)
    [~, baseName] = fileparts(files(k).name);

    % Full input path without extension 
    fname_image = fullfile(srcFolder, baseName);

    % White reference
    fname_whiteReference = fullfile(pwd, '20011415_smooth.xlsx');

    % Running processing
    [im_reflectance, WL] = hsf_rad2reflectance(fname_image, fname_whiteReference);

   
% Saving  all outputs into  folder, with unique names
    outFile = fullfile(dstFolder, baseName + "_reflectance.mat");
    save(outFile, 'im_reflectance', 'WL', '-v7.3');
end
