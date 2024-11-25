clear, clc, close all

points_per_virus = [1; 2; 5; 10; 20]

for l = 1:size(points_per_virus);

for f = 1:10;

N_DIM   = 2; % 2D
D       = 0.05; % 0.1 - 0.01 µm^2/s --> 0.032 for HLADR
D_slow  = D/points_per_virus(l,1); % 0.015 for H18+HLADR
dT      = 0.03; % s,% Time step between acquisition; fast acquisition!
time_steps = 1000;   

%% Generate random point in a square 

%clear, clc, close all

% Set the number of points
N_PARTICLES = 100;

% Generate random x and y coordinates within the 1x1 µm square
x = rand(N_PARTICLES, 1); % 100 random numbers for x coordinates
y = rand(N_PARTICLES, 1); % 100 random numbers for y coordinates

% Combine x and y coordinates into a single matrix
starting_points = [x, y];

% Display the points
scatter(x, y, 'filled');
xlabel('x (µm)');
ylabel('y (µm)');
title('100 Random Points in a 1x1 µm Square');
axis([0 1 0 1]); % Set axis limits to match the 1x1 µm square
axis square
grid on;

% Generate Viruses

N_VIRUSES = 10;
radius = 0.05;
radii = radius * ones(N_VIRUSES, 1);

% Generate random x and y coordinates within the 1x1 µm square
x = (1-radius) * rand(N_VIRUSES, 1); % 100 random numbers for x coordinates
y = (1-radius) * rand(N_VIRUSES, 1); % 100 random numbers for y coordinates

% Combine x and y coordinates into a single matrix
viruses = [x, y];

viscircles(viruses, radii)

%% Step-by-step

close all, clc

k       = sqrt(2 * D * dT);
k_vir   = sqrt(2 * D_slow * dT);

tracks       = cell(N_PARTICLES, 1);
step = 1;

for i = 1 : N_PARTICLES

    % Initial position
    trajectory = [];
    X0 = starting_points(i,:);
    position = X0;
    trajectory(1, :) = position;

    % Simulate the random walk
    step = 1;
    
    while step<time_steps;
        
    step = size(trajectory,1)+1;    
    distances = sqrt(sum((viruses - position).^2, 2));

    if sum(distances<radius)==0 % if its not in one of the circles
           
    % Store the position
    displacement = k * randn(1, N_DIM);  % Generate random displacements
    position = position + displacement;     % Update the position
    %trajectory(step, :) = position;
    
    else
    displacement = k_vir * randn(1, N_DIM);
    position = position + displacement;    
    %trajectory(step, :) = position;
    end
       
    if or(position(1,1)<0, position(1,1)>1)==1
        
        trans_vector = trajectory(step-1, 1:2) - position;
        new_position = trajectory(step-1, 1:2) + trans_vector;
        position = new_position;    
        trajectory(step, :) = position;
        
    elseif or(position(1,2)<0, position(1,2)>1)==1
        
        trans_vector = trajectory(step-1, 1:2) - position;
        new_position = trajectory(step-1, 1:2) + trans_vector;
        position = new_position;    
        trajectory(step, :) = position;
    else 
    
       trajectory(step, :) = position;
    end
       
    end
    
    clear position X0  
    
    % Store the trajectory
    time = [];    
    time = (0 : size(trajectory,1)-1)' * dT;
    tracks{i} = [time trajectory];
 
end

clear i X dX time X0


ID = 1; 
plot(tracks{ID, 1}(:,2),tracks{ID, 1}(:,3), 'k-', 'LineWidth', 2); hold on
viscircles(viruses, radii)
axis equal;
% xlim([-1, 1]);
% ylim([-1, 1]);
xlabel('X');
ylabel('Y');
title('Random Points within a Circle');
grid on;
hold off;

%% Plot all trajectories

close all

particles_per_virus = []; distances = []

figure
for ID = 1:N_PARTICLES;
    
plot(tracks{ID, 1}(:,2),tracks{ID, 1}(:,3), 'b-', 'LineWidth', 2); hold on

end
viscircles(viruses, radii);
axis square;

xlabel('X');
ylabel('Y');
title('Random Points within a Circle');
grid on;
hold off;

end_points = [];

figure
for ID = 1:N_PARTICLES;
    
scatter(tracks{ID, 1}(end,2),tracks{ID, 1}(end,3),'filled'); hold on

end_points(ID,:) = [tracks{ID, 1}(end,2) tracks{ID, 1}(end,3)];

end
viscircles(viruses, radii);
axis square
title(['Median = ' num2str(median(particles_per_virus))])

% Count number of points in each circle

for i = 1:N_VIRUSES;
   
  distances = sqrt(sum((end_points - viruses(i,:)).^2, 2));
  particles_per_virus(i,1) = sum(sum(distances<radius));

end

figure

hist(pdist(end_points),50)
title(['Median = ' num2str(median(pdist(end_points)))])

points_per_virus(l,f+1) = mean(particles_per_virus);

%close all

end

end



%% Make animated movie of one moving particle - Step 1 Check ID

ID = 2; 

plot3(x_circle, y_circle,z_circle_base, 'b-', 'LineWidth', 2); hold on
plot3(x_circle, y_circle,z_circle_top, 'b-', 'LineWidth', 2); hold on
plot3(x_points(ID), y_points(ID),z_points(ID), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); hold on
plot3(tracks{ID, 1}(:,2),tracks{ID, 1}(:,3),tracks{ID, 1}(:,4), 'g-', 'LineWidth', 2); hold on

axis equal;
xlim([-radius, radius]);
ylim([-radius, radius]);
zlim([heigth, -heigth/2]);
xlabel('X');
ylabel('Y');
title(['Track length ' num2str(length(tracks{ID}))]);
grid on;
hold off;

%% Step 2 animate

figure('Position', [600 600 400 300])
set(gcf, 'color', 'w')

% Initialize view angle
az = 90;
el = 5;  % Set a fixed elevation angle

h = plot3(tracks{ID, 1}(:,2),tracks{ID, 1}(:,3),tracks{ID, 1}(:,4), 'g-', 'LineWidth', 2); hold on;
m = plot3(tracks{ID, 1}(:,2),tracks{ID, 1}(:,3),tracks{ID, 1}(:,4), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Set axis limits and aspect ratio
axis([-radius radius -radius radius heigth -heigth/2]);
axis vis3d;  % Prevents MATLAB from changing the aspect ratio when rotating
axis equal;
xlabel('X [µm]');
ylabel('Y [µm]');
zlabel('Z [µm]');
title('Random Points within a Circle');

% Initialize handle for the timestamp text
time_text_handle = text(-1.4, 1.3, 1.4, '', 'FontSize', 12, 'Color', 'red');

% Initialize GIF
filename = 'diffusion_animation.gif';  % Output GIF file name
frame = getframe(gcf);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);


for i = 1:60:size(tracks{ID},1);
    
plot3(x_circle, y_circle,z_circle_base, 'b-', 'LineWidth', 2); 
plot3(x_circle, y_circle,z_circle_top, 'b-', 'LineWidth', 2); 
plot3(x_points(ID), y_points(ID),z_points(ID), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 

% Rotate the view
    az = az + 1;  % Increment azimuth angle for rotation
    view(az, el);  % Update view angle
   
    
set(h, 'XData', tracks{ID, 1}(1:i,2), 'YData',tracks{ID, 1}(1:i,3), 'ZData', tracks{ID, 1}(1:i,4));
    drawnow;  % Force MATLAB to update the figure window
    pause(0.01);  % Pause for a short time (in seconds)

    
set(m, 'XData', tracks{ID, 1}(i,2), 'YData',tracks{ID, 1}(i,3), 'ZData', tracks{ID, 1}(i,4));
    drawnow;  % Force MATLAB to update the figure window
    pause(0.01);  % Pause for a short time (in seconds) 
    
    % Update timestamp
    time_text = sprintf('Elapsed time: %d min', i);title(time_text, 'FontSize', 14);
%     if ishandle(time_text_handle)
%         delete(time_text_handle);  % Delete the previous timestamp
%     end
%     time_text_handle = text(-1.4, 1.4, 1.4, time_text, 'FontSize', 12, 'Color', 'red');
    

 % Capture the plot as an image and write to GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    
    pause(0.01);  % Pause for a short time (in seconds)    
    
end


