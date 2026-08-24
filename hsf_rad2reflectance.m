function [im_reflectance, WL] = hsf_rad2reflectance(fname_image, fname_whiteReference)
% Validating input paths before reading the image and reference file.
    if nargin < 2 || isempty(fname_image) || isempty(fname_whiteReference)
        error('hsf_rad2reflectance:InputError', ...
            'Both fname_image and fname_whiteReference must be provided.');
    end

    imgFile = [fname_image, '.img'];

    if ~exist(imgFile, 'file')
        error('hsf_rad2reflectance:FileNotFound', 'Image file not found: %s', imgFile);
    end

    if ~exist(fname_whiteReference, 'file')
        error('hsf_rad2reflectance:FileNotFound', 'White reference file not found: %s', fname_whiteReference);
    end

    % Selecting the ROI and collecting the white-reference rows.
    [im_cropped, ~, im_white, ~, ~] = hsf_crop(imgFile);

    % Reading the reference spectrum from the spreadsheet.
    WL_reference = xlsread(fname_whiteReference, 'Sheet1', 'A3:A2253');
    R_reference  = xlsread(fname_whiteReference, 'Sheet1', 'B3:B2253');

    % Reading the wavelength grid from the ENVI header.
    info = hsi_readinfo(fname_image);
    WL = info.lambda(1:end-4)';

    % Interpolating the white-reference reflectance onto the image wavelength grid.
    Rwhite = interp1(WL_reference, R_reference, WL, 'linear');
    im_reflectance = hs_raw2ref(im_cropped(:,:,1:end-4), im_white(:,:,1:end-4), Rwhite);

    % Displaying the mean reflectance and white-reference radiance plots.
    figure;
    subplot(2,2,[1,3]);
    imagesc(mean(im_reflectance, 3));
    axis equal tight;
    caxis([0 100]);
    colorbar;
    title('Mean Reflectance');

    subplot(2,2,2);
    plot(WL, mean(unfold(im_reflectance))');
    title('\color{magenta}Mean image spectra');
    ylabel('Reflectance');
    xlabel('Wavelength');

    subplot(2,2,4);
    plot(squeeze(mean(im_white,2))');
    title('\color{green}White reference radiance');
    ylabel('Radiance');
    xlabel('N Column');
end

function im_reflectance = hs_raw2ref(im_cropped, im_white, Rwhite)
% Normalizing the sample radiance by the white reference and applying the reference curve.
    mean_im_white = mean(im_white, 1);
    imR = im_cropped ./ repmat(mean_im_white, size(im_cropped, 1), 1, 1);
    imR_uf = unfold(imR);

    % Flattening each pixel spectrum for element-wise multiplication by the white spectrum.
    Rwhite_imRuf = repmat(Rwhite, size(imR, 1) * size(imR, 2), 1);
    im_reflectance_uf = imR_uf .* Rwhite_imRuf;
    im_reflectance = reshape(im_reflectance_uf, size(imR,1), size(imR,2), size(imR,3));
end

function Xuf = unfold(X)
% Flattening the first two dimensions while keeping the spectral dimension intact.
    [x, y, z] = size(X);
    Xuf = reshape(X, x * y, z);
end

function [im_cropped, mask, im_white, rgb_cropped, coord] = hsf_crop(fname, varargin)
% Selecting the crop window and extracting the white-reference rows from the same columns.
    [pathstr, name, ~] = fileparts(fname);
    header_fname = fullfile(pathstr, name);
    info = hsi_readinfo(header_fname);

    if nargin == 1
        X = multibandread(fname, [info.lines info.samples info.bands], ...
            info.datatype, info.offset, info.interleave, 'ieee-le', ...
            {'band', 'Direct', info.defaultbands});
        rgb = X;

        imshow(mean(rgb,3), [prctile(rgb(:),10) prctile(rgb(:),90)]);
        title('Select the area to crop and DBLCLK inside');
        [rgb_cropped, mask, coord] = hs_crop(rgb);
    else
        if islogical(varargin{1})
            mask = varargin{1};
            [coord(:,2), coord(:,1)] = find(mask > 0);
            rgb_cropped = [];
        else
            rgb = varargin{1};
            imagesc(uint16(rgb));
            [rgb_cropped, mask, coord] = hs_crop(rgb);
        end
    end

    CY1 = min(coord(:,2));
    CY2 = max(coord(:,2));
    CX1 = min(coord(:,1));
    CX2 = max(coord(:,1));

    % Reading only the selected spatial region from the hyperspectral cube.
    im_cropped = multibandread(fname, [info.lines info.samples info.bands], info.datatype, ...
        info.offset, info.interleave, 'ieee-le', ...
        {'row', 'range', [CY1, CY2]}, ...
        {'column', 'range', [CX1, CX2]});

    if nargin > 2
        lines_for_white = varargin{2};
    else
        lines_for_white = [10, 30];
    end

    im_white = multibandread(fname, [info.lines info.samples info.bands], info.datatype, ...
        info.offset, info.interleave, 'ieee-le', ...
        {'row', 'range', [lines_for_white(1), lines_for_white(end)]}, ...
        {'column', 'range', [CX1, CX2]});

    rectangle('Position', [CX1, lines_for_white(1), CX2 - CX1, lines_for_white(2) - lines_for_white(1)], 'EdgeColor', 'g');
    title('\color{magenta}Cropped region; \color{green}white reference');
end

function info = hsi_readinfo(fname)
% Reading the ENVI header, mapping data types, and extracting wavelength metadata.
    info = [];

    d = {'uint8', 'int16', 'int32', 'single', 'double', 'uint16', 'uint32', 'int64', 'uint64'};
    ds = [1 2 4 4 8 2 4 8 8 8];

    if ~ischar(fname)
        error('fname should be a char string');
    end

    rfid = fopen(strcat(fname, '.hdr'), 'r');
    if rfid == -1
        error('error while opening %s.hdr', fname);
    end

    t = [];
    while true
        tline = fgetl(rfid);
        if ~ischar(tline)
            break;
        end

        [first, second] = strtok(tline, '=');
        first(first == ' ') = [];

        switch lower(first)
            case 'wavelength'
                while isempty(find(second == '}'))
                    second = [second fgetl(rfid)];
                end
                [~, s] = strtok(second);
                s(find(s == '{')) = ' ';
                s(find(s == '}')) = ' ';
                info.lambda = strread(s, '%f', 'delimiter', ',');

            case 'interleave'
                [~, s] = strtok(second);
                s(find(s == ' ')) = [];
                info.interleave = lower(s);

            case 'headeroffset'
                [~, s] = strtok(second);
                info.offset = str2num(s);

            case 'samples'
                [~, s] = strtok(second);
                info.samples = str2num(s);

            case 'lines'
                [~, s] = strtok(second);
                info.lines = str2num(s);

            case 'bands'
                [~, s] = strtok(second);
                info.bands = str2num(s);

            case 'datatype'
                [~, s] = strtok(second);
                t = str2num(s);

            case 'defaultbands'
                while isempty(find(second == '}'))
                    second = [second fgetl(rfid)];
                end
                [~, s] = strtok(second);
                s(find(s == '{')) = ' ';
                s(find(s == '}')) = ' ';
                info.defaultbands = strread(s, '%f', 'delimiter', ',');

                switch t
                    case 1
                        t = d{1};
                        sz = ds(1);
                    case 2
                        t = d{2};
                        sz = ds(2);
                    case 3
                        t = d{3};
                        sz = ds(3);
                    case 4
                        t = d{4};
                        sz = ds(4);
                    case 5
                        t = d{5};
                        sz = ds(5);
                    case 12
                        t = d{6};
                        sz = ds(6);
                    case 13
                        t = d{7};
                        sz = ds(7);
                    case 14
                        t = d{8};
                        sz = ds(8);
                    case 15
                        t = d{9};
                        sz = ds(9);
                    otherwise
                        error('Unknown image data type');
                end
                info.datatype = t;
                info.datasize = sz;
        end
    end

    fclose(rfid);
end

function [imc, mask, coord] = hs_crop(im)
    h = imrect;
    pos0 = wait(h);
    setColor(h, 'm');

    pos = round(pos0);
    x1 = pos(2);
    x2 = pos(2) + pos(4);
    y1 = pos(1);
    y2 = pos(1) + pos(3);

    imc = im(x1:x2, y1:y2, :);
    mask = zeros(size(im, 1), size(im, 2));
    mask(x1:x2, y1:y2) = 1;

    coord = zeros((x2 - x1 + 1) * (y2 - y1 + 1), 2);
    k = 0;
    for i = x1:x2
        for j = y1:y2
            k = k + 1;
            coord(k,:) = [j, i];
        end
    end
end