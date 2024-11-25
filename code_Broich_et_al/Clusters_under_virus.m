clear, clc, close all

%% 1. Load data

folder = 'S:\GENERAL\Primary_Lukas\Experiments\lbr_Imaging\lbr_0017\Analysiert\Data\X31\006\Analysis';
cd(folder)
filename = 'All_Trajectories';
output_folder = 'S:\GENERAL\Primary_Lukas\Experiments\lbr_Imaging\lbr_0017\Analysiert\Data\X31\006\Analysis';

WF_folder = 'S:\GENERAL\Primary_Lukas\Experiments\lbr_Imaging\lbr_0017\Analysiert\Data\X31\006';
WF_name = 'C2-2022_09_29_lbr0017_EGFR-mEos_MDCK_X31_WF_006.tif';

locs_input = dlmread([filename '.csv'],',',1,0);

% Find the respective Columns

file = fopen([filename '.csv']); % csv for TS
header = fgetl(file);
header = regexp(header, ',', 'split');

xCol                = strmatch('"x [nm]"',header);
yCol                = strmatch('"y [nm]"',header);
zCol                = strmatch('"z [nm]"',header);
photonsCol          = strmatch('"intensity [photons]"',header);
framesCol           = strmatch('"frame"',header);
SigmaCol            = 6 %strmatch('"sigma_x [nm]"',header);

fprintf('\n -- Data loaded --\n')


xCol = 1; % in pxl
yCol = 2;
framesCol = 4;
photonsCol = 3;
SigmaCol = 8;

locs_input(:,1) = locs_input(:,1)*160;
locs_input(:,2) = locs_input(:,2)*160;

%% Load WF image

cd(WF_folder)
WF = imread(WF_name);
% imshow(WF)
WF2 = flipdim(imrotate(WF,-90),2);

imshow(flipdim(imrotate(WF,-90),2));hold on;

% imagesc(flipdim(imrotate(WF,-90),2));hold on;

%% Filtering and overlay with WF
close all
minFrame            = 1000;
MinPhotons          = 100;
maxSig              = 30;
        
filter              = [];
filter              = find(locs_input(:,photonsCol) > MinPhotons ...       % Photon filter
                         & locs_input(:,framesCol) > minFrame ...          % Frame filter
                         & locs_input(:,SigmaCol) < maxSig);               % Sigma filter
                                              
locs            = locs_input(filter,1:end);



figure('units','normalized','outerposition',[0 0 1 1])
imagesc(WF2,[0 20000]);hold on;
colormap gray
% imagesc(imadjust(imgaussfilt(rendered,1)));hold on;
% scatter(locs(:,xCol)./160,locs(:,yCol)./160,'.','red');
% imshow(imgaussfilt(rendered,1),[0.5 100]);
% colormap hot
scatter(locs(:,xCol)./160,locs(:,yCol)./160,5,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.05)

axis square
box on

%% 2. Peform DBSCAN clustering (with manual ROI)
close all

clear Rg

cd('S:\GENERAL\Primary_Christian\Stella\analysis_scripts\cluster\cluster\functions');

pxlsize = 160;

k       = 30; 
Eps     = 20;

[subset, FiC, CpA, LpC, T] = DBSCAN_with_ROI_coloc(WF2,locs, xCol, yCol, pxlsize, k, Eps);

for i = 1:max(subset(:,end));
    
    target  = find(subset(:,end)==i);
    cluster = subset(target,1:end);
    
    Rg(i,1)  = sqrt(sum(var([cluster(:,xCol), cluster(:,yCol)],1,1))); % Rg 2D
    Ecc(i,1) = sqrt(max(eig(cov(cluster(:,xCol:yCol))))/min(eig(cov(cluster(:,xCol:yCol))))); % Ecc XY

end

ROI = 8;

cd(output_folder)
save([filename '_Rg_ROI_' num2str(ROI) '.mat'],'Rg');
save([filename '_Ecc_ROI_' num2str(ROI) '.mat'],'Ecc');
save([filename '_T_ROI_' num2str(ROI) '.mat'],'T');

median(Rg)

%% 

figure('Position',[700 500 300 100])


Result = [Rg Ecc];
T = table(Result);

uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


%% Correlate 

% Control
       
rendered = hist3([locs(:,yCol),locs(:,xCol)],[128 128]);
figure
imshow(imgaussfilt(rendered,1),[0.5 50]);hold on;
scatter(locs(:,xCol)./160,locs(:,yCol)./160,'.','red');
%% 


CorrX = []; CorrY = [];

for x = 2:129;
    
    for y = 2:129;

% target = find((locs(:,xCol))>(x-1)*160 & (locs(:,xCol))< x*160 & (locs(:,yCol))>(y-1)*160 & (locs(:,yCol)) < y*160);

target = find((locs(:,xCol))>(x-1)*160 & (locs(:,xCol))< (x+0)*160 & (locs(:,yCol))>(y-1)*160 & (locs(:,yCol)) < (y+0)*160);


%CorrX = horzcat(CorrX, WF(y-1,x-1));
CorrX = horzcat(CorrX, rendered(y-1,x-1));
CorrY = horzcat(CorrY, size(target,1));

    end
end

figure
% scatter(CorrX, CorrY,'.')
% CorrY = cast(CorrY,'uint8');

% CorrX = im2double(CorrX);
scatter(CorrX/max(CorrX), CorrY/max(CorrY),'.')
R = corrcoef(CorrX/max(CorrX), CorrY/max(CorrY));

title(['My ' num2str(min(R))])
min(R)

%% Test

close all

x = 10; y = 60;

WF3 = imgaussfilt(rendered,1);
WF3(y-1,x-1) = 256;

figure
% imshow(WF3);hold on;
imshow(imgaussfilt(WF3,1),[0.5 50]);hold on;

target = find((locs(:,xCol))>(x-2)*160 & (locs(:,xCol))< (x+1)*160 & (locs(:,yCol))>(y-2)*160 & (locs(:,yCol)) < (y+1)*160);
scatter(locs(target,xCol)./160,locs(target,yCol)./160,'.','red');
title(['My ' num2str(size(target,1))])

