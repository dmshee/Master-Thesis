function [im_reflectance, WL]=hsf_rad2reflectance(fname_image,fname_whiteReference)
% function to calculate reflectance from HYSPEX RADIANCE FILES.
% Calculates reflectance by normalizing a radiance image: im_cropped, by a radiance image of white reference: im_white;  and correct by the white
% reference reflectance spectra.Ana Herrero-Langreo,23/09/2021, HMM, UCD.
% 
% INPUT:
% filename of excel file with white reference
% fname_whiteReference='20011415_smooth.xlsx' % 100% 
%
% filename of hyperspectral image WITHOUT EXTENSION
% fname_image='BACTERIAALL_10OD_AL_REP1REP2_SWIR_384_SN3168_8000us_2020-12-17T122604_raw_rad'
%
% OUTPUT:
% im_reflectance> corrected reflectance image of the cropped area
% WL> wavebands
%
% 
% [im_reflectance, WL]=hsf_rad2reflectance(fname_image,fname_whiteReference)
%
%% load radiance image and ask if we want to select a region of interest
close all,
[im_cropped,mask,im_white,rgb_cropped,coord]=hsf_crop([fname_image,'.img']);

%% load white reference spectra to be used to correct the image
WL_reference=xlsread(fname_whiteReference,'Sheet1','A3:A2253');
R_reference=xlsread(fname_whiteReference,'Sheet1','B3:B2253');

%% calculate reflectance (Reflectance sample =( Radiance sample ./ Radiance white spectralon placed with the sample) .* reference spectra of the white spectralon from xls file
% interpolate white reference reflectance
info=hsi_readinfo(fname_image); % info from header file
WL=info.lambda(1:end-4)';% I exclude the last 4 wavelenths, they are problematic for SWIR

Rwhite = interp1(WL_reference,R_reference,WL,'linear');% interpolated reference spectra of the spectralon (from xls file)

[im_reflectance]=hs_raw2ref(im_cropped(:,:,1:end-4),im_white(:,:,1:end-4),Rwhite);

figure,
subplot(2,2,[1,3]),imagesc(mean(im_reflectance,3)), axis equal tight, caxis([0 100]),colorbar
title('Mean Reflectance')
subplot(2,2,2), plot(WL,mean(unfold(im_reflectance))'), title('\color{magenta}Mean image Spectra'), ylabel('Reflectance'),xlabel('Wavelength');
subplot(2,2,4),plot(squeeze(mean(im_white,2))'),title('\color{green}White Reference Radiance'),ylabel('Radiance'),xlabel('N Column');
end

function [im_reflectance]=hs_raw2ref(im_cropped,im_white,Rwhite)
% function to calculate reflectance f.rom files already LOADED in the
% workspace
% It calculates reflectance by normalizing a radiance image: im_cropped, by a radiance image of white reference: im_white;  and correct by the white
% reference reflectance spectra.Ana Herrero-Langreo,18/02/2020, HMM, UCD.
% 
%[im_reflectance]=hs_raw2ref(im_cropped,im_white,Rwhite)
 
 mean_im_white=mean(im_white,1);
imR=im_cropped./repmat(mean_im_white,size(im_cropped,1),1,1);
imR_uf=unfold(imR);
Rwhite_imRuf=repmat(Rwhite,size(imR,1)*size(imR,2),1);
im_reflectance_uf=imR_uf.*Rwhite_imRuf;
im_reflectance=reshape(im_reflectance_uf,size(imR,1),size(imR,2),size(imR,3));
end
 
function[Xuf]=unfold(X)
%[Xuf]=unfold(X)
[x,y,z]=size(X);
Xuf=reshape(X,x*y,z);
end

function [im_cropped,mask,im_white,rgb_cropped,coord]=hsf_crop(fname,varargin)
% function to crop hyperspectral image from file. Loads only the frames of
% the image specified, opens an imagesc of those planes, asks the user to
% select the roi of the image and loads all wavebands for the roi selected.
% %as well as corresponding white reference at the beginning of the image.
%Ana Herrero-Langreo, UCD,
% 2020, HMM project, written for 16uint precision .Hyspex files. Imshow
% does not work on uint32 on matlab 2018
%
% INPUT:
% fname: filename of hyperspectral image with extension i.e:
% fname='EC_REP3_IND_VNIR_1800_SN00863_12500us_2020-01-31T124352_raw_rad.img'
% OPTIONAL INPUTS:
% mask (binary mask to cut the image)
% rgb (rgb or plane), opens rgb or plane image to manually select the area
% to crop
%
% OUTPUTS:
% im_cropped: cropped hyperspectral image
% mask 
% im_white: 50 first pushbroom lines of the image at the columns corresponding with the cropped area, normally corresponding
% to white reference
% coord: row and column coordinates of the selected area
%
% [im_cropped,mask,im_white,rgb_cropped,coord]=hsf_crop(fname)% to
% generate RGB image and manually select area to crop.
%
% [im_cropped,mask,im_white,rgb_cropped,coord]=hsf_crop(fname,rgb)% to
% manually select area to crop from RGB image.
%
% [im_cropped,mask,im_white,rgb_cropped,coord]=hsf_crop(fname,mask)% to
% crop image using a mask (mask has to be in "logical" format)
%
% [im_cropped,mask,im_white,rgb_cropped,coord]=hsf_crop(fname,[],[10,30])% to
%  selects lines 10 to 30 as white reference,generate RGB image and manually select area to crop. 
%
% it assumes .hdr file exists!


% extract info from header file
[pathstr,name,ext] = fileparts(fname);
header_fname = fullfile(pathstr,name);
info=hsi_readinfo(header_fname);
%% rgb or hsi plane used to cut the image: 
if nargin==1 % if only the filename is provided, obtains a rgb image and plots it to manually crop the hsi
    X = multibandread(fname,[info.lines info.samples info.bands],info.datatype,info.offset,info.interleave,'ieee-le', {'band','Direct',info.defaultbands});
    %rgb=uint16(X);
    rgb = X;
    
    %imshow(rgb,[]);
    imshow(mean(rgb,3),[prctile(rgb(:),10) prctile(rgb(:),90)]);
    title('Select the area to crop and DBLCLK inside')
    [rgb_cropped,mask,coord]=hs_crop(rgb); 
else
    if islogical (varargin{1}) % if mask is provided, to automatically crop hsi around the mask
    mask=varargin{1};
    [coord(:,2),coord(:,1)]=find(mask>0);
    rgb_cropped=[];
    else
    rgb=varargin{1}; % if rgb (or hyperspectral plane(s)) is provided, to manually crop hsi
    imagesc(uint16(rgb))
    [rgb_cropped,mask,coord]=hs_crop(rgb);    
    end
end

%% extract coordinates to cut the image

CY1=min(coord(:,2));
CY2=max(coord(:,2));
CX1=min(coord(:,1));
CX2=max(coord(:,1));
% figure,imagesc(rgb(CY1:CY2,CX1:CX2,1))

% read hs image only within the coordinates
im_cropped = multibandread(fname,[info.lines info.samples info.bands],info.datatype,...
    info.offset,info.interleave,'ieee-le',{'row','range',[CY1,CY2]},...
    {'column','range',[CX1,CX2]});

% extract image of corresponding white reference at the same columns as the
% cut image.
if nargin>2
   lines_for_white=varargin{2};
else
lines_for_white=[10,30]
end
im_white= multibandread(fname,[info.lines info.samples info.bands],info.datatype,...
    info.offset,info.interleave,'ieee-le',{'row','range',[lines_for_white(1),lines_for_white(end)]},...
    {'column','range',[CX1,CX2]}); %  first lines of the image,specified by "lines_for_white" white reference
rectangle('Position',[CX1,lines_for_white(1),CX2-CX1,lines_for_white(2)-lines_for_white(1)],'EdgeColor','g')
title('\color{magenta}Cropped region; \color{green}white reference')
end
%%
%%
function info=hsi_readinfo( fname )
% info=hsi_readinfo( ficname );
%
% read the header file nomfic.hdr 
%
% return informations into a structure :
%  .samples    : number of samples 
%  .lines      : number of lines 
%  .bands      : number of bands
%  .offset     : offset of binary file
%  .datatype   : type of data (bit8, ..., uint64)
%  .datasize   : size in bytes of data
%  .interleave : type of interleave
%  .lambda     : wavelengths of the bands
% .defaultbands: bands for reconstructing rgb image
% JM ROGER, UMR ITAP, Cemagref, 2008

info=[];

% Parameters initialization
elements={'samples' 'lines' 'bands' 'datatype' 'headeroffset' 'interleave','defaultbands'};
d={'uint8' 'int16' 'int32' 'single' 'double' 'uint16' 'uint32' 'int64' 'uint64'};
ds=[1 2 4 4 8 2 4 8 8 8];

% Check user input
if ~ischar(fname)
    error('fname should be a char string');
end


% Open ENVI header file to retreive s, l, b & d variables
rfid = fopen(strcat(fname,'.hdr'),'r');

if rfid==-1
    error(sprintf('error while opening %s.hdr', fname));
    return;
end;

% Read ENVI image header file and get p(1) : nb samples,
% p(2) : nb lines, p(3) : nb bands and t : data type

while 1
    tline = fgetl(rfid);
    
    if ~ischar(tline), break, end
    
    [first,second]=strtok(tline,'=');
    first(first==' ')=[];
    
    switch lower(first)
        case 'wavelength'
            
            while isempty(find(second=='}'))
                second = [second fgetl(rfid)];
            end;
            [f,s]=strtok(second);
            s(find(s=='{'))=' ';
            s(find(s=='}'))=' ';
            info.lambda = strread( s, '%f', 'delimiter',',' );
        case 'interleave'
            [f,s]=strtok(second);
            s(find(s==' '))=[];
            info.interleave = lower(s);
        case 'headeroffset'
            [f,s]=strtok(second);
            info.offset = str2num(s);
        case 'samples'
            [f,s]=strtok(second);
            info.samples = str2num(s);
        case 'lines'
            [f,s]=strtok(second);
            info.lines = str2num(s);
        case 'bands'
            [f,s]=strtok(second);
            info.bands = str2num(s);
        case 'datatype'
            [f,s]=strtok(second);
            t=str2num(s);
        case 'defaultbands'
             while isempty(find(second=='}'))
                second = [second fgetl(rfid)];
            end;
            [f,s]=strtok(second);
            s(find(s=='{'))=' ';
            s(find(s=='}'))=' ';
            info.defaultbands = strread( s, '%f', 'delimiter',',' );
             
            switch t
                case 1
                    t=d{1};
                    sz=ds(1);
                case 2
                    t=d{2};
                    sz=ds(2);
                case 3
                    t=d{3};
                    sz=ds(3);
                case 4
                    t=d{4};
                    sz=ds(4);
                case 5
                    t=d{5};
                    sz=ds(5);
                case 12
                    t=d{6};
                    sz=ds(6);
                case 13
                    t=d{7};
                    sz=ds(7);
                case 14
                    t=d{8};
                    sz=ds(8);
                case 15
                    t=d{9};
                    sz=ds(9);
                otherwise
                    error('Unknown image data type');
            end
            info.datatype=t;
            info.datasize=sz;
    end
end
fclose(rfid);
end
%%
function [imc,mask,coord]=hs_crop(im)
%function to crop a hyperspectral image from pinpoint in the screen,needs to have a plane or an hsi
%figure open.
% written by Ana Herrero-langreo (ana.herrero-langreo@ucd.ie), for preliminar exp. CAC2018, proj. HMM, 28-02-2018 
%
% INPUT
% im: hs image
%
% OUTPUT
%imc: cropped hs image
%mask: from the original image, 1: cropped image, 0: discarded bacground
%
% imagesc(im(:,:,25));
% hs_show (im);
% [imc,mask]=hs_crop(im);


% [CX,CY] = ginput(2);
% pos=round([CX,CY]);
% p11=pos(1,2); p12=pos(2,2); p21=pos(1,1);p22=pos(2,1);
h=imrect;
pos0=wait(h);

setColor(h,'m');
pos=round(pos0);
x1=pos(2);x2=(pos(2)+pos(4));y1=pos(1);y2=(pos(1)+pos(3));

imc=im(x1:x2,y1:y2,:);
mask=zeros(size(im,1),size(im,2));
mask(x1:x2,y1:y2)=1;

 
 coord = zeros((x2-x1+1)*(y2-y1+1),2);
 k=0;
 for i=x1:x2
     for j=y1:y2
         k=k+1;
         coord(k,:)=[j,i];
     end;
 end;

end