close all 
clear 
clc
%% Parameters 
% Sporangiophore parameters 
r_cylinder = 1;         % Sporangiophore radius 
H0 = 10*r_cylinder;                 % Initial length of the cylinder 
H_final = 1.8*H0;         % Final length of the cylinder 

% Sporangium parametrs 
r_cap = 1.5*r_cylinder;   % radius of the sporangium 

% Rotation Parameter
rot_ang = 2*pi; 

% Arrow distribution on the sporangium
theta_dist = [80 65 40 25 10]*pi/180; 
L_dist = [1/6 2/6 3/6 4/6 5/6]; 

% Animation parameters 
N_frames = 50;               % Number of frames for the animation
num_arrows = 8; 

% Frame Number (for single frame plot) 
Frame_number = 25; 
%% % Set up the figure
fig1 = figure(1);
set(gcf,'Color','white')
set(gcf,"Units","inches")
hold on;
axis equal;
view(3);
axis([-r_cap r_cap -r_cap r_cap 0 (H_final + r_cylinder)]);
axis off

%% Cylinder and the Spherical cap 
% Create cylinder data
thetavec = linspace(0, 2*pi, 50);
[Theta, Z] = meshgrid(thetavec, linspace(0, 1, 50)); % Z scaled from 0 to 1 for easy height adjustment
X_cylinder = r_cylinder * cos(Theta);
Y_cylinder = r_cylinder * sin(Theta);

% Create spherical cap data (height equal to radius)
cap_height = r_cylinder; 
cap_angle = pi - acos(1 - (cap_height / r_cap)); % Calculate angle based on height
[theta_cap, phi_cap] = meshgrid(linspace(0, 2*pi, 50), linspace(0, cap_angle, 25));
X_cap = r_cap * sin(phi_cap) .* cos(theta_cap);
Y_cap = r_cap * sin(phi_cap) .* sin(theta_cap);
Z_cap = r_cap * cos(phi_cap);


% Plot initial objects 
cylinder_surf = surf(X_cylinder, Y_cylinder, Z * H0, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
cap_surf = surf(X_cap, Y_cap, Z_cap + H0,'FaceColor', [0, 0, 0], 'EdgeColor', 'none'); 

%% Create initial arrows on the cylinder 
% Define arrows on the surface of the cylinder
arrow_length = 1*(2*pi*r_cylinder/num_arrows);           % Length of each rotation arrow
arrow_theta = linspace(0, 2*pi, num_arrows + 1); % Angles for arrows, evenly spaced
arrow_theta(end) = [];        % Remove the last point to avoid overlap at 0 and 2*pi
arrows = gobjects(length(theta_dist), num_arrows); 


for mm = 1 : length(theta_dist)
    % Each following Loop makes the initial arrows at one height along the
    % sporagiophore 
    for kk = 1 : num_arrows
        % Arrows at Z = 2 (initial middle position, red color)
        arrow_X = r_cylinder * cos(arrow_theta(kk));
        arrow_Y = r_cylinder * sin(arrow_theta(kk));
        arrow_Z = H0 * L_dist(mm); % Middle position
        arrows(mm, kk) = quiver3(arrow_X, arrow_Y, arrow_Z, ...
                               arrow_length * cos(arrow_theta(kk) + pi/2), ...
                               arrow_length * sin(arrow_theta(kk) + pi/2), ...
                               arrow_length * sin(theta_dist(mm)), 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
    end
end



%% Animation loop
for i = 1 : N_frames-2
    % Calculate interpolated height for the cylinder
    current_height = H0 + (H_final - H0) * (i / N_frames);
    Z_cylinder_scaled = Z * current_height;

    % Calculate rotation angle (full 360 degrees over N_frames)
    angle = rot_ang * (i / N_frames); % Rotation angle in radians

    % Apply rotation to the X and Y coordinates of the cylinder and cap
    X_cylinder_rotated = X_cylinder * cos(angle) - Y_cylinder * sin(angle);
    Y_cylinder_rotated = X_cylinder * sin(angle) + Y_cylinder * cos(angle);

    X_cap_rotated = X_cap * cos(angle) - Y_cap * sin(angle);
    Y_cap_rotated = X_cap * sin(angle) + Y_cap * cos(angle);

    % Update cylinder and cap data with the rotated and scaled height
    set(cylinder_surf, 'XData', X_cylinder_rotated, 'YData', Y_cylinder_rotated, 'ZData', Z_cylinder_scaled);
    set(cap_surf, 'XData', X_cap_rotated, 'YData', Y_cap_rotated, 'ZData', Z_cap + current_height);

    % Update arrow positions and directions to indicate rotation and elongation
       for nn = 1 : length(theta_dist)
        for kk = 1:num_arrows
            % Update arrows at Z = 2 (initial middle position, red)
            arrow_X_rotated = r_cylinder * cos(arrow_theta(kk) + angle);
            arrow_Y_rotated = r_cylinder * sin(arrow_theta(kk) + angle);
            set(arrows(nn, kk), 'XData', arrow_X_rotated, 'YData', arrow_Y_rotated, ...
                             'ZData', current_height * L_dist(nn)); % Middle position based on current height
            set(arrows(nn, kk), 'UData', arrow_length * cos(arrow_theta(kk) + angle + pi/2), ...
                             'VData', arrow_length * sin(arrow_theta(kk) + angle + pi/2), ...
                             'WData', arrow_length * sin(theta_dist(nn))); % Increased elongation component for visual effect
    
        end
       end
    % Capture the current frame for the GIF
    drawnow; % Ensure the figure is updated
    frame = getframe(gcf); % Capture the current figure

end


%%  Set up the second figure for initial configuration only 
% Same line of codes are copied before the loop section 
fig2 = figure(2);
set(gcf,'Color','white')
set(gcf,"Units","inches")
hold on;
axis equal;
view(3);
axis([-r_cap r_cap -r_cap r_cap 0 (H_final + r_cylinder)]);
axis off

%% Cylinder and the Spherical cap 
% Create cylinder data
thetavec = linspace(0, 2*pi, 50);
[Theta, Z] = meshgrid(thetavec, linspace(0, 1, 50)); % Z scaled from 0 to 1 for easy height adjustment
X_cylinder = r_cylinder * cos(Theta);
Y_cylinder = r_cylinder * sin(Theta);

% Create spherical cap data (height equal to radius)
cap_height = r_cylinder; 
cap_angle = pi - acos(1 - (cap_height / r_cap)); % Calculate angle based on height
[theta_cap, phi_cap] = meshgrid(linspace(0, 2*pi, 50), linspace(0, cap_angle, 25));
X_cap = r_cap * sin(phi_cap) .* cos(theta_cap);
Y_cap = r_cap * sin(phi_cap) .* sin(theta_cap);
Z_cap = r_cap * cos(phi_cap);


% Plot initial objects 
cylinder_surf = surf(X_cylinder, Y_cylinder, Z * H0, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
cap_surf = surf(X_cap, Y_cap, Z_cap + H0,'FaceColor', [0, 0, 0], 'EdgeColor', 'none'); 

%% Create initial arrows on the cylinder 
% Define arrows on the surface of the cylinder
arrow_length = 1*(2*pi*r_cylinder/num_arrows);           % Length of each rotation arrow
arrow_theta = linspace(0, 2*pi, num_arrows + 1); % Angles for arrows, evenly spaced
arrow_theta(end) = [];        % Remove the last point to avoid overlap at 0 and 2*pi
arrows = gobjects(length(theta_dist), num_arrows); 


for mm = 1 : length(theta_dist)
    % Each following Loop makes the initial arrows at one height along the
    % sporagiophore 
    for kk = 1 : num_arrows
        % Arrows at Z = 2 (initial middle position, red color)
        arrow_X = r_cylinder * cos(arrow_theta(kk));
        arrow_Y = r_cylinder * sin(arrow_theta(kk));
        arrow_Z = H0 * L_dist(mm); % Middle position
        arrows(mm, kk) = quiver3(arrow_X, arrow_Y, arrow_Z, ...
                               arrow_length * cos(arrow_theta(kk) + pi/2), ...
                               arrow_length * sin(arrow_theta(kk) + pi/2), ...
                               arrow_length * sin(theta_dist(mm)), 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
    end
end

%% Animation loop
for i = 1 : Frame_number
    % Calculate interpolated height for the cylinder
    current_height = H0 + (H_final - H0) * (i / N_frames);
    Z_cylinder_scaled = Z * current_height;

    % Calculate rotation angle (full 360 degrees over N_frames)
    angle = rot_ang * (i / N_frames); % Rotation angle in radians

    % Apply rotation to the X and Y coordinates of the cylinder and cap
    X_cylinder_rotated = X_cylinder * cos(angle) - Y_cylinder * sin(angle);
    Y_cylinder_rotated = X_cylinder * sin(angle) + Y_cylinder * cos(angle);

    X_cap_rotated = X_cap * cos(angle) - Y_cap * sin(angle);
    Y_cap_rotated = X_cap * sin(angle) + Y_cap * cos(angle);

    % Update cylinder and cap data with the rotated and scaled height
    set(cylinder_surf, 'XData', X_cylinder_rotated, 'YData', Y_cylinder_rotated, 'ZData', Z_cylinder_scaled);
    set(cap_surf, 'XData', X_cap_rotated, 'YData', Y_cap_rotated, 'ZData', Z_cap + current_height);

    % Update arrow positions and directions to indicate rotation and elongation
       for nn = 1 : length(theta_dist)
        for kk = 1:num_arrows
            % Update arrows at Z = 2 (initial middle position, red)
            arrow_X_rotated = r_cylinder * cos(arrow_theta(kk) + angle);
            arrow_Y_rotated = r_cylinder * sin(arrow_theta(kk) + angle);
            set(arrows(nn, kk), 'XData', arrow_X_rotated, 'YData', arrow_Y_rotated, ...
                             'ZData', current_height * L_dist(nn)); % Middle position based on current height
            set(arrows(nn, kk), 'UData', arrow_length * cos(arrow_theta(kk) + angle + pi/2), ...
                             'VData', arrow_length * sin(arrow_theta(kk) + angle + pi/2), ...
                             'WData', arrow_length * sin(theta_dist(nn))); % Increased elongation component for visual effect
    
        end
       end
    % Capture the current frame for the GIF
    drawnow; % Ensure the figure is updated
    frame = getframe(gcf); % Capture the current figure

end
