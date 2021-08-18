%% Load the localization file 

clear, clc, close all

Locpath1         = ['S:\GENERAL\Primary_Christian\01_data\Nikon_TIRF\decode'];
locName1         = ['emitter_newModel'];

cd(Locpath1);
locs_Ch1        =  readtable([locName1 '.csv']);

pxlsize = 70; % 70 nm for pco.edge 4.2

fprintf('\n -- Data Loaded --\n');

%% RCC Drift Correct

% Generate Coords variable as input for RCC
% Remove NaN and Inf

pxlsize = 70; % 70 nm for pco.edge 4.2

clc
cd([fitting_dist '/RCC_drift_correct']);

coords(:,1) = locs(:,xCol)/pxlsize;
coords(:,2) = locs(:,yCol)/pxlsize;
coords(:,3) = locs(:,framesCol);

% Remove NaN and Inf
temp = coords;
clear coords
coords = temp( ~any( isnan( temp(:,1) ) | isinf( temp(:,1) ), 2 ),: );

fprintf('\n -- Ready for DC -- \n')


%% Select region to correct

pxlsize = 1;

heigth  = round((max(coords(:,2))-min(coords(:,2)))/pxlsize);
width   = round((max(coords(:,1))-min(coords(:,1)))/pxlsize);
im      = hist3([coords(:,1),coords(:,2)],[width heigth]); % heigth x width

% Select rectangles

rect = []; 

figure('Position',[100 200 600 600])
f = imagesc(imrotate(im,90),[0 2e3]);
colormap('parula'); colorbar;

rect = getrect;

fprintf('\n -- ROI selected --\n')

xmin = min(coords(:,1))+ rect(1,1)*pxlsize;
ymin = max(coords(:,2))- rect(1,2)*pxlsize - (rect(1,4)*pxlsize) ;
xmax = xmin + 256;
ymax = ymin + 256;

target      = find(coords(:,1)>xmin & coords(:,1)<xmax & coords(:,2)>ymin & coords(:,2)<ymax);
coords_ROI  = coords(target,1:end);

% Show cropped region

heigth  = round((max(coords_ROI(:,2))-min(coords_ROI(:,2)))/pxlsize);
width   = round((max(coords_ROI(:,1))-min(coords_ROI(:,1)))/pxlsize);
im      = hist3([coords_ROI(:,1),coords_ROI(:,2)],[width heigth]); % heigth x width

figure('Position',[100 200 600 600])
imagesc(imrotate(im,90),[0 2e3]);
colormap('parula'); colorbar;

coords_ROI(:,1) = coords_ROI(:,1)-min(coords_ROI(:,1));
coords_ROI(:,2) = coords_ROI(:,2)-min(coords_ROI(:,2));


%% Drift correct
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

segpara     = max(locs(:,framesCol))/5; 
imsize      = 256;
pixelsize 	= 106;
binsize     = 30;

[coordscorr, finaldrift] = RCC(coords_ROI, segpara, imsize, pixelsize, binsize, 0.2);
% [coordscorr, finaldrift] = DCC(coords, segpara, imsize, pixelsize, binsize);
                          
clc

display(['DC finished in ' round(num2str(toc/60)) ' min']);

%% Plot Drift curves

figure('Position',[100 100 900 400])
subplot(2,1,1)
plot(finaldrift(:,1))
title('x Drift')
subplot(2,1,2)
plot(finaldrift(:,2))
title('y Drift')

%% Apply correction to coords

deltaXY = [];
deltaXY(:,1) = coords(:,3); % frames
deltaXY(:,2) = finaldrift(deltaXY(:,1),1); % x drift
deltaXY(:,3) = finaldrift(deltaXY(:,1),2); % y drift

coordsDC      = [];
coordsDC(:,1) = coords(:,1)-deltaXY(:,2);
coordsDC(:,2) = coords(:,2)-deltaXY(:,3);

display('Coord corrected');
   
% Generate DC variable

locs_DC = locs;
locs_DC(:,xCol) = coordsDC(:,1)*pixelsize; % x nm
locs_DC(:,yCol) = coordsDC(:,2)*pixelsize; % y nm

close all

%% Generate output for ThunderSTORM

tic;
cd(Locpath1);;

locsTS = [];
locsTS(:,1) = locs_Ch1.frame_ix;          % frames
locsTS(:,2) = locs_Ch1.x*pxlsize;               % x nm
locsTS(:,3) = locs_Ch1.y*pxlsize;               % y nm
locsTS(:,4) = locs_Ch1.z*pxlsize;               % z nm
locsTS(:,5) = locs_Ch1.phot;         % photons
locsTS(:,6) = locs_Ch1.x_sig*pxlsize;               % LL

NameCorrected = [locName1 '_decode.csv'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[['"frame","x [nm]","y [nm]","z [nm]","intensity [photons]", "sigma [nm]" '] ' \n']);
dlmwrite(NameCorrected,locsTS,'-append');
fclose('all');

clc;

clc
display(['File saved in ' num2str(toc) ' min']);




