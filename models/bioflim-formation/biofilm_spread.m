function [u, h_bio] = biofilm_spread(r, s, u_init)
% Hybrid spread model (spread of surface density and thickness)
% 1D along-catheter length model with coupled surface biomass u (CFU/cm^2)
% and biofilm thickness (cm)
% - u evolves by reaction (logistic), lateral diffusion, conversion to volume
% - h grows by converted biomass (and detaches)

%% Parameters
% Inputs
D = 0.01;         % diffusion coefficient (cm^2/hr)
k_conv = 0.1;    % conversion rate (1/hr) from surface biomass -> accumulated volume
rho = 1e11;      % volumetric packing density (CFU / cm^3)
u_max = 1e7;     % maximum density/capacity (CFU/cm^2)

% Space parameters --> x in [x0, x1]
h_space = 0.005;    % space discretization step (dx)
x0 = 0;             % left boundary
x1 = s;             % right boundary (s = catheter length (cm))
x = x0:h_space:x1;  % spatial grid
J = length(x);      % number of spatial points (size of space grid)

%% Initial Conditions
u = zeros(J, 1);       % 2D biofilm grid: thickness x length
newu = zeros(J, 1);    % next time step
h_bio = zeros(J, 1);   % thickness (cm) -- initially zero

% Initialize biofilm density with gaussian bump
ind = x < s/10;
u(ind) = u_max * exp(-x(ind).^2 / 0.2);

% 2nd way to initialize bump
% sigma = 0.5;  % bump width in cm
% u = u_init * exp(-(x.^2)/(2*sigma^2))';  % Gaussian bump

% Initial thickness consistent with initial u
h_bio = u / rho;    % (cm)

%% Laplacian Matrix
% Laplacian discretization with 'no flux' Neumann boundary conditions
e = ones(J, 1);
L = spdiags([e -2*e e], -1:1, J, J);
% Enforce "no-flux" B.C. by modifying first and last rows (ghost method)
L(1, :) = 0;
L(J, :) = 0;

%% Time Parameters
t0 = 0;
tfinal = 35;    % units: hours

% Time Discretization
dt = min(0.5, 0.9 * (h_space^2 / 2/D));  % CFL condition
Nt = ceil((tfinal - t0) / dt);     % number of time steps

% Preallocate time-series storage
t_values = zeros(Nt, 1);
avg_u = zeros(Nt, 1);
avg_h = zeros(Nt, 1);

% Snapshot storage for overlays
snapshot_interval = 5;
num_snapshots = ceil(tfinal / snapshot_interval) + 2;
snapshots_u = zeros(num_snapshots, J);
snapshots_h = zeros(num_snapshots, J);
snap_times = zeros(num_snapshots, 1);
snap_count = 0;

% Plot Setup
% biofilm (u) profile along catheter length (x)
figure(1);
clf;
subplot(2, 1, 1); % overlayed spatial profiles (u)
hold on;
title('Spatial profiles (overlay)');
xlabel('Position along catheter (cm)');
ylabel('Biofilm density (CFU/cm^2)');

subplot(2, 1, 2); % overlayed thickness profile
hold on;
title('Spatial profile: thickness');
xlabel('Position along catheter (cm)');
ylabel('Thickness (\mum)')

% uncomment the 2 lines below to see the ic before running the simulation
%disp('')
%pause

%% Main Loop (explicit Euler)
t = t0;
step_index = 1;
tic                 % start a stopwatch timer
while t < tfinal
    % Reaction + diffusion + conversion (explicit)
    % Note: u in CFU/cm^2, diffusion uses spatial units in metres
    Lu = (L * u);     % discrete laplacian w/ interior points
    newu = u + dt * (r .* u.*(1 - u ./ u_max) + (D / h_space^2) * Lu);

    % Enforce Neumann BC (no flux): set boundary equal to neighbour
    newu(1) = newu(2);        % enforce no-flux at left (u(x0) = u(x0+h))
    newu(J) = newu(J - 1);    % enforce no-flux at right (u(x1) = u(x1-h))

    % Thickness update (cm): accumulation from conversion
    newh = h_bio + dt * (k_conv ./ rho) .* u;

    % Update
    u = newu;
    u = max(0, min(u, u_max));   % ensure u stays in [0, u_max)
    h_bio = newh;

    % Record time series
    t_values(step_index) = t;
    avg_u(step_index) = mean(u);
    avg_h(step_index) = mean(h_bio);

    % Save snapshot @ each interval
    if (abs(mod(t, snapshot_interval)) < 1e-8) || (t == 0) || (t+dt > tfinal)
        snap_count = snap_count + 1;
        snapshots_u(snap_count, :) = u';
        snapshots_h(snap_count, :) = h_bio';
        snap_times(snap_count) = t;
    end

    % Live overlay plotting
    if (step_index == 1) || (abs(mod(t, 0.1)) < 1e-8)
        subplot(2, 1, 1);
        plot(x, u, 'LineWidth', 1);
        axis([x0 x1 0 u_max*1.1]);
        drawnow;

        subplot(2, 1, 2);
        plot(x, h_bio*1e4, 'LineWidth', 1); % convert thickness to microns
        axis([x0 x1 0 max(h_bio)*1e4*1.1+1e-6]);    % fix max y bound to 100 microns
        drawnow;
    end

    % Advance time
    t = t + dt;
    step_index = step_index + 1;

    fprintf("t = %.5f\n", t);

end
toc     % stops stopwatch & display elapsed time since tic was called

% Trim unused snapshot rows
snapshots_u = snapshots_u(1:snap_count, :);
snapshots_h = snapshots_h(1:snap_count, :);
snap_times = snap_times(1:snap_count);

%% Visualization - Overlay all Snapshots
figure;
subplot(2, 1, 1);
hold on;
for k = 1:snap_count
    plot(x, snapshots_u(k, :), 'DisplayName',sprintf('t = %.1f', snap_times(k)));
end
xlabel('Position along catheter (cm)');
ylabel('Biofilm density (CFU/cm^2)');
axis([x0 x1 0 u_max*1.1]);
legend('show');

subplot(2, 1, 2);
hold on;
for k = 1:snap_count
    plot(x, snapshots_h(k, :) * 1e4, 'DisplayName', sprintf('t = %.1f', snap_times(k)));
end

xlabel('Position along catheter (cm)');
ylabel('Thickness (\mum)');
legend('show');

% Plot time series of averages
figure;
yyaxis left
plot(t_values(1:step_index - 1), avg_u(1:step_index - 1), 'LineWidth', 1.5);
ylabel('Average u (CFU/cm^2)');
ylim([0 u_max*1.1]);

yyaxis right
plot(t_values(1:step_index - 1), avg_h(1:step_index - 1) * 1e4, 'LineWidth', 1.5);
ylabel('Average thickness (\mum)');
xlabel('Time');
end