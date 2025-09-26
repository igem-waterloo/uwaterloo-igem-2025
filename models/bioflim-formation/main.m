% Main file
close all
clear
clc

%% Variables
catheter_length = 40;     % catheter length (cm)
D = 5e-6;                 % diffusion coefficient of bacteria (cm^2/s)

%% Biofilm Attachment
[cells, time, avg] = biofilm_attachment(catheter_length, D);
    % outputs are in the units: million cells/cm^2, seconds
formation_rate_time = time;    % time to 30% coverage (sec)

% Convert formation rate and initial film density to CFU/cm^2 to input into
% the biofilm spread model
formation_rate_cells = cells * 1e6 / 3600;  % CFU/cm^2/hr
avg_cell_density = avg * 1e6;        % CFU/cm^2

%% Biofilm Spread
% Spread is measured in CFU/cm^2 over along the catheter length over time
[spread, thickness] = biofilm_spread(formation_rate_cells, catheter_length, avg_cell_density);

%% Biofilm Detachment
growth_velocity = sqrt(cells * 1e6) * 1e-3 / 2 / 3600; % cm/hr
biofilm_detachment_main(growth_velocity);
