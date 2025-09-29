function [] = biofilm_detachment_main(growth_velocity)
%% Constants 
D = 0.1;         % diffusion coefficient (cm^2/hr)
gif_filename = 'biofilm_detachment.gif';

%% Sptial Parameters
h_space = 0.05;    % space discretization step (dx)
s = 40; % catheter length (cm)
x0 = 0; % left boundary
x1 = s; % right boundary
x = x0:h_space:x1; % sptial grid
J = length(x);
u_max = 1e7;

%% Time Parameters
t0 = 0;
tfinal = 35; % hours
t = t0;

dt = min(0.5, 0.9 * (h_space^2 / 2/D)); % hours
Nt = ceil((tfinal - t0) / dt);     % number of time steps

%% Plot Setup
% Biofilm thickness (L) along catheter length (x)
figure(4);
clf;
title('Spatial profiles (overlay)');
xlabel('Position along catheter (cm)');
ylabel('Biofilm Thickness (m)');
lines_plotted = 0;
step_index = 1;

%% Main Loop
L = zeros(J,1);
density = zeros(J, 1); 
ind = x < s/10;
density(ind) = u_max * exp(-x(ind).^2 / 0.2);
while t < tfinal
    [density, thickness] = biofilm_spread_for_detachment(0.5, 0.1, density);
    density(isnan(density)) = 0;

    L = biofilm_detachment(density, t, L, dt, growth_velocity);
    L(isnan(L)) = 0;

     % Live overlay plotting
    
    if ((step_index == 1) || (abs(mod(t, 0.1)) < 1e-8))
        plot(x, L, 'LineWidth', 1);
        axis([x0 x1 0 1e-5]);
        title('Spatial profiles (overlay)');
        xlabel('Position along catheter (cm)');
        ylabel('Biofilm Thickness (m)');
         % Capture the frame for the GIF
        frame = getframe(gcf);  % Capture current figure as a frame
        [A, map] = rgb2ind(frame.cdata, 256);  % Convert to indexed image format
        
        % Write the first frame and subsequent frames to the GIF
        if t == t0  % First frame (creates GIF file)
            imwrite(A, map, gif_filename, 'gif', 'LoopCount', inf, DelayTime = 0.7);
        else
            imwrite(A, map, gif_filename, 'gif', WriteMode='append', DelayTime = 0.4);  % Append to existing GIF
        end
        drawnow;
        lines_plotted = lines_plotted + 1;
          
    end
    t = t + dt;
    step_index = step_index + 1;
    fprintf("t = %.5f\n", t); 
end