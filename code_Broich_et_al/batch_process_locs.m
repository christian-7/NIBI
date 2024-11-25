clear, clc, close all

Input_folder             = 'S:\GENERAL\Primary_Christian\01_data\Nikon_NSTORM\2024-05-14_Ninj2_Ella\locResults\locResults_mEos';
Output_folder            = 'S:\GENERAL\Primary_Christian\01_data\Nikon_NSTORM\2024-05-14_Ninj2_Ella\analysis';

cd(Input_folder);
image_files = dir(sprintf('*.h5'));

%% Analyze the localization file

% Open the loc file

for file_ID = 1:size(image_files,1);
    
clear locs_xyz   locs_frames locs_photons locs_xyz_sig locs_prob locs_bg coords
 
cd(image_files(file_ID).folder);
h5disp([image_files(file_ID).name]);

locs_xyz        = transpose(h5read(image_files(file_ID).name,'/data/xyz'));
locs_frames     = transpose(h5read(image_files(file_ID).name,'/data/frame_ix'));
locs_photons    = transpose(h5read(image_files(file_ID).name,'/data/phot'));
locs_xyz_sig    = transpose(h5read(image_files(file_ID).name,'/data/xyz_sig'));
locs_prob       = transpose(h5read(image_files(file_ID).name,'/data/prob'));
locs_bg         = transpose(h5read(image_files(file_ID).name,'/data/bg'));

pxlsize = 70; % 70 nm for pco.edge4.2, 160 nm for Andor iXon

fitting_dist = 'S:\GENERAL\Primary_Christian\02_Image_analysis\spline_fitting_distribution';

fprintf('\n -- Data Loaded --\n');

cd ..

%% DC

% RCC Drift Correct

% Generate Coords variable as input for RCC
% Remove NaN and Inf

% pxlsize = 160; % 70 nm for pco.edge 4.2

clc
cd([fitting_dist '/RCC_drift_correct']);

coords(:,1) = locs_xyz(:,1);
coords(:,2) = locs_xyz(:,2);
coords(:,3) = locs_frames+1;

% Remove NaN and Inf
temp = coords;
clear coords
coords = temp( ~any( isnan( temp(:,1) ) | isinf( temp(:,1) ), 2 ),: );

fprintf('\n -- Ready for DC -- \n')

% Drift correct
tic
close all
% Input:    coords:             localization coordinates [x y frame], 
%           segpara:            segmentation parameters (time wimdow, frame)
%           imsize:             image size (pixel)
%           pixelsize:          camera pixel size (nm)
%           binsize:            spatial bin pixel size (nm)
%           rmax:               error threshold for re-calculate the drift (pixel)
% Output:   coordscorr:         localization coordinates after correction [xc yc] 
%           finaldrift:         drift curve (save A and b matrix for other analysis)

%finaldrift = RCC_TS(filepath, 1000, 256, 160, 30, 0.2);

segpara     = max(coords(:,3))/5; 
imsize      = 256;
pixelsize 	= 70;
binsize     = 30;

[coordscorr, finaldrift] = RCC(coords, segpara, imsize, pixelsize, binsize, 0.2);
% [coordscorr, finaldrift] = DCC(coords, segpara, imsize, pixelsize, binsize);
                          
clc

display(['DC finished in ' round(num2str(toc/60)) ' min']);

% Plot Drift curves

figure('Position',[100 100 900 400])
subplot(2,1,1)
plot(finaldrift(:,1))
title('x Drift')
subplot(2,1,2)
plot(finaldrift(:,2))
title('y Drift')

% Apply correction to coords
   
%% Generate DC variable

locs_DC = locs_xyz;
locs_DC(:,1) = coordscorr(:,1)*pxlsize; % x nm
locs_DC(:,2) = coordscorr(:,2)*pxlsize; % y nm

close all

display('DC finished');

% Filter localizations

locs_DC(:,4)    = locs_frames;
locs_DC(:,5)    = locs_photons;
locs_DC(:,6)    = locs_prob;
locs_DC(:,7:8)  = locs_xyz_sig(:,1:2)*pxlsize;
locs_DC(:,9)    = locs_xyz_sig(:,3);

minFrame            = 100;
MinPhotons          = 100; % S/P: 500/100
MaxPhotons          = 10000;
MinProb             = 0.8;
z_min               = -1000;
z_max               = 1000; 
maxSig              = [50 50 10000];
        
filter              = [];
filter              = find(locs_DC(:,5) > MinPhotons & locs_DC(:,5) < MaxPhotons ...                        % Photon filter
                         & locs_DC(:,6) > MinProb ...                                                       % Prob filter
                         & locs_DC(:,4) > minFrame ...                                                      % Frame filter
                         & locs_DC(:,3) > z_min & locs_DC(:,3) < z_max ...                                  % Z range filter
                         & locs_DC(:,7) < maxSig(1) & locs_DC(:,8) < maxSig(2) & locs_DC(:,9) < maxSig(3)); % Sigma filter
                                              
locsFilt            = locs_DC(filter,1:end);

clc
display(['Localizations filtered (' num2str(length(locsFilt)/length(locs_DC)) ' left)']);

f = msgbox(['Localizations filtered (' num2str(length(locsFilt)/length(locs_DC)) ' left)']);
set(f, 'position', [500 500 200 50]); %makes box bigger
th = findall(f, 'Type', 'Text');                   %get handle to text within msgbox
th.FontSize = 12;


% Render localizations

%locsFilt = locs_DC; % xCol = 1; yCol = 2;

pxlsize_render = 10; % nm

heigth = round((max(locsFilt(:,2)) - min((locsFilt(:,2))))/pxlsize_render);
width  = round((max(locsFilt(:,1)) - min((locsFilt(:,1))))/pxlsize_render);
        
rendered = hist3([locsFilt(:,2),locsFilt(:,1)],[heigth width]);

% figure 
% imagesc(imadjust(imgaussfilt(rendered,1)));
% imshow(imgaussfilt(rendered,1),[0.5 100]);
% colormap hot


%% Save Tiff file

cd(Output_folder);

I32 = [];
I32 = uint32(rendered);

name = [strrep(image_files(file_ID).name,'.h5','') '_10nm_pxl.tiff'];

t = Tiff(name,'w');
tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()


%% Generate output for ThunderSTORM

tic;
cd(Output_folder);

locsTS = [];
locsTS(:,1) = locsFilt(:,4);               % frames
locsTS(:,2) = locsFilt(:,1);               % x nm
locsTS(:,3) = locsFilt(:,2);               % y nm
locsTS(:,4) = locsFilt(:,3);               % z nm
locsTS(:,5) = locsFilt(:,5);               % photons
locsTS(:,6) = locsFilt(:,7);               % sigma x
locsTS(:,7) = locsFilt(:,8);               % sigma y
locsTS(:,8) = locsFilt(:,9);               % sigma z

NameCorrected = [strrep(image_files(file_ID).name,'.h5','') '_decode.csv'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[['"frame","x [nm]","y [nm]","z [nm]","intensity [photons]", "sigma_x [nm]", "sigma_y [nm]", "sigma_z [nm]" '] ' \n']);
dlmwrite(NameCorrected,locsTS,'-append');
fclose('all');

clc;

clc
display(['File saved in ' num2str(toc) ' min']);

end

