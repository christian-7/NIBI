%% Perform clustering of SMLM data
% 
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% 1. Load sample data
% 2. Calculate Ripley's K statistics
% 3. Perform Ripley's K clustering

% Ripleys K
% This approach compares the measured distribution of single-molecule localizations to a simulated random distribution, and provides information whether clustering occurs. 
% The maximum of the H function reflects the average size of clusters. 
% The amplitude of the H-function is a measure of the degree of clustering. 
% The variance of the amplitude reflects the variance of the spatial organization of proteins within a cluster.


% Dependencies: 
% 
% DBSCAN_with_ROI
% manualROI
% calculateRipleysK
% calculateLKfunction
% DBSCAN

clear, clc, close all

%% 1. Load data

folder = 'S:\GENERAL\Primary_Lukas\Experiments\P_Paper Project Immobilization\Data for Paper\Tracking\PR8\006\Analysis';
cd(folder)
filename = 'A549-EGFRmEos_PR8-NiR_006_filtered_localizations';
output_folder = 'S:\GENERAL\Primary_Lukas\Experiments\P_Paper Project Immobilization\Data for Paper\Tracking\PR8\clustering_for_revision';

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

%% Filtering

minFrame            = 500;
MinPhotons          = 110;
maxSig              = 0.3;
        
filter              = [];
filter              = find(locs_input(:,photonsCol) > MinPhotons ...       % Photon filter
                         & locs_input(:,framesCol) > minFrame ...          % Frame filter
                         & locs_input(:,SigmaCol) < maxSig);               % Sigma filter
                                              
locs            = locs_input(filter,1:end);

pxlsize_render = 20; % nm

heigth = round((max(locs(:,yCol)) - min((locs(:,yCol))))/pxlsize_render);
width  = round((max(locs(:,xCol)) - min((locs(:,xCol))))/pxlsize_render);
        
rendered = hist3([locs(:,yCol),locs(:,xCol)],[heigth width]);

figure 
imagesc(imadjust(imgaussfilt(rendered,1)));
%imshow(imgaussfilt(rendered,1),[0.5 100]);
colormap hot
axis square

%% 2. Peform DBSCAN clustering (with manual ROI)
close all
cd('S:\GENERAL\Primary_Christian\Stella\analysis_scripts\cluster\cluster\functions');

pxlsize = 200;

k       = 30; %20 number of points
Eps     = 20; %40 radius

%[subset, FiC, CpA, LpC, T] = DBSCAN_with_ROI(locs, xCol, yCol, pxlsize, k, Eps);
[subset,locs_ROI, FiC, CpA, LpC, T] = DBSCAN_with_ROI_buildin(locs, xCol, yCol, pxlsize, k, Eps);

for i = 1:max(subset(:,end));
    
    target  = find(subset(:,end)==i);
    cluster = subset(target,1:end);
    
    Rg(i,1)  = sqrt(sum(var([cluster(:,xCol), cluster(:,yCol)],1,1))); % Rg 2D
    Ecc(i,1) = sqrt(max(eig(cov(cluster(:,xCol:yCol))))/min(eig(cov(cluster(:,xCol:yCol))))); % Ecc XY

end

ROI = 9;

cd(output_folder)
save([filename '_Rg_ROI_' num2str(ROI) '.mat'],'Rg');
save([filename '_Ecc_ROI_' num2str(ROI) '.mat'],'Ecc');
save([filename '_T_ROI_' num2str(ROI) '.mat'],'T');

% Save the localizations

save([filename '_Clusters_' num2str(ROI) '.mat'],'subset');

%% 3. Calculate Ripley's K statistics (with manual ROI)
% Generate loc file for LAMA *_Lama.txt
 
close all
cd('S:\GENERAL\Primary_Christian\Stella\analysis_scripts\cluster\cluster\functions');
 
[locs_ROI] = manualROI(locs,xCol,yCol,200);
 
ROI     = [min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) ... 
           min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))];
 
maxk    = 1000;    % maximum distance
stepk   = 10;     % step size
 
[K, K_rand] = calculateRipleysK(locs_ROI, xCol, yCol, ROI, maxk, stepk);
 
L = sqrt(K(:,2)/pi);
H = L - K(:,1);
M = K(:,2)./(pi*(K(:,1).^2));
 
L_rand = sqrt(K_rand(:,2)/pi);
H_rand = L_rand - K_rand(:,1);
M_rand = K_rand(:,2)./(pi*(K_rand(:,1).^2));
 
figure('Position',[100 400 500 500])
 
subplot(2,2,1)
plot(K(:,1),K(:,2),'-b');hold on
plot(K_rand(:,1),K_rand(:,2),'-r');hold on
ylabel('K(r)'); box on; axis square; xlabel('distance [nm]');
 
subplot(2,2,2)
plot(K(:,1),L,'b');hold on
plot(K_rand(:,1),L_rand,'r');hold on
ylabel('L(r)');box on; axis square; xlabel('distance [nm]');
 
subplot(2,2,3)
plot(K(:,1),H,'b');hold on
plot(K_rand(:,1),H_rand,'r');hold on
ylabel('L(r)-r');box on; axis square; xlabel('distance [nm]');
 
subplot(2,2,4)
plot(K(:,1),M,'b');hold on
plot(K_rand(:,1),M_rand,'r');hold on
ylabel('K(r)/(\pi*r^2)'); box on; axis square; xlabel('distance [nm]');


%%  Save L curve

Var1(:,1) = K(:,1);
Var1(:,2) = H;
Var1(:,3) = K_rand(:,1);
Var1(:,4) = H_rand;

ROI = 1;

cd(output_folder)
save([filename '_Lrr_ROI_' num2str(ROI) '.mat'],'Var1');

%% 4. Calculate Ripley's K and L statistics (with manual ROI)

% Code from: Shivanandan A, et al (2015) PLOS ONE 10(3): e0118767. https://doi.org/10.1371/journal.pone.0118767
% https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118767

close all

radius = 5:5:500;

[locs_ROI] = manualROI(locs,xCol,yCol,200);

ROI     = [min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) ... 
           min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))];

ROI_cropped = find(locs_ROI(:,xCol)>ROI(1)+max(radius) & locs_ROI(:,xCol)<ROI(2)-max(radius) & ...
                   locs_ROI(:,yCol)>ROI(3)+max(radius) & locs_ROI(:,yCol)<ROI(4)-max(radius));


[TrueLrr, TrueK] = calculateLKfunction(locs_ROI(:,1:2), radius);

[max_value,index] = max(TrueLrr)

figure('Position',[100 400 600 300])

subplot(1,2,1)
plot(radius,TrueK,'-b');hold on
ylabel('K(r)'); box on; axis square; xlabel('distance [nm]');

subplot(1,2,2)
plot(radius,TrueLrr,'b');hold on
title(['Maximum at ' num2str(max(radius(index))) ' nm'])
ylabel('L(r)-r');box on; axis square; xlabel('distance [nm]');

%% Save L curve

Var1(:,1) = radius;
Var1(:,3) = TrueLrr;

ROI = 1;

cd(output_folder)
save([filename '_Lrr_ROI_' num2str(ROI) '.mat'],'Var1');


