function [formation_rate_cells, formation_rate_time, avg_cell_density] = biofilm_attachment(catheter_length, D_bacteria)
    % e. faecalis biofilm attachment model
    % returns formation rate in million cells/cm²/s and time to 30%
    % coverage (critical point for clinical considerations
    
    %% Parameters
    %catheter_length = 40;     % catheter length (cm)
    catheter_diameter = 0.5;  % inner diameter (cm)
    T = 37;                  % body temperature (°C)
    urine_flow_rate = 1.5;   % mL/min
    
    cross_area = pi * (catheter_diameter/2)^2;
    velocity = (urine_flow_rate/60) / cross_area; % cm/s
    
    % model parameters
    ka_protein = 0.08;       % protein adsorption rate (/s)
    kd_protein = 0.005;      % protein desorption rate (/s)
    C_protein = 2.0;         % protein concentration (mg/mL)
    %D_bacteria = 5e-6;       % bacterial diffusion coefficient (cm²/s)
    n_bulk = 1e7;           % bacterial concentration (cells/mL)
    ka_attach_base = 0.12;   % base attachment rate (/s)
    
    %biofilm density parameters
    max_cell_density = 5e8;  % maximum cells per cm² in mature biofilm
    cell_height = 1e-4;      % average cell height (cm) for E. faecalis
    
    %% setting up arrays
    dx = 5.0;
    x = 0:dx:catheter_length;
    nx = length(x);
    
    t_max = 3600;  % 1 hour
    dt = 30;       % 30 second steps
    time = 0:dt:t_max;
    nt = length(time);
    
    theta_protein = zeros(nx, nt);
    theta_biofilm = zeros(nx, nt);
    cell_density = zeros(nx, nt);  % cells per cm²
    
    %% main simulation
    for t_idx = 1:nt
        for x_idx = 1:nx
            current_x = x(x_idx);
            
            % conditioning film formation with langmuir adsorption model
            if t_idx == 1
                theta_protein(x_idx, t_idx) = 0;
            else
                dtheta_dt = ka_protein * C_protein * (1 - theta_protein(x_idx, t_idx-1)) - ...
                           kd_protein * theta_protein(x_idx, t_idx-1);
                theta_protein(x_idx, t_idx) = theta_protein(x_idx, t_idx-1) + dtheta_dt * dt;
                theta_protein(x_idx, t_idx) = max(0, min(1, theta_protein(x_idx, t_idx)));
            end
            
            % bacterial transport, advection-diffusion eqn
            residence_time = current_x / velocity;
            mass_transfer_coeff = D_bacteria / (catheter_diameter/4);
            transport_efficiency = 1 - exp(-mass_transfer_coeff * residence_time);
            n_surface = n_bulk * transport_efficiency;
            
            % biofilm final attachment, smulchowski model
            attachment_efficiency = 0.2 + 0.6 * theta_protein(x_idx, t_idx);
            position_factor = exp(-current_x / 20);
            ka_attach = ka_attach_base * attachment_efficiency * position_factor;
            
            if t_idx == 1
                theta_biofilm(x_idx, t_idx) = 0;
                cell_density(x_idx, t_idx) = 0;
            else
                attachment_flux = ka_attach * (n_surface/n_bulk) * ...
                                (1 - theta_biofilm(x_idx, t_idx-1));
                theta_biofilm(x_idx, t_idx) = theta_biofilm(x_idx, t_idx-1) + ...
                                            attachment_flux * dt;
                theta_biofilm(x_idx, t_idx) = max(0, min(1, theta_biofilm(x_idx, t_idx)));
                
                cell_density(x_idx, t_idx) = theta_biofilm(x_idx, t_idx) * max_cell_density;
            end
        end
    end
    
    %% formation rates calculation

    avg_biofilm = mean(theta_biofilm, 1);
    idx_30 = find(avg_biofilm >= 0.3, 1, 'first');
    
    if isempty(idx_30)
        formation_rate_time = inf;
    else
        formation_rate_time = time(idx_30);
    end
   
    avg_cell_density = mean(cell_density, 1); % average cells/cm² over catheter length
    
    % rate of change in cell density (million cells/cm²/s)
    if length(time) > 1
        density_rate = diff(avg_cell_density) ./ diff(time); % cells/cm²/s
        formation_rate_cells = mean(density_rate) / 1e6; % convert to million cells/cm²/s
        
        % rate during active attachment phase (first 30 minutes)
        active_phase_idx = time <= 1800; % first 30 minutes
        if sum(active_phase_idx) > 1
            active_density = avg_cell_density(active_phase_idx);
            active_time = time(active_phase_idx);
            active_rate = diff(active_density) ./ diff(active_time);
            formation_rate_cells = mean(active_rate) / 1e6; % million cells/cm²/s
        end
    else
        formation_rate_cells = 0;
    end
    
    % Update average cell density for output
    avg_cell_density = avg_cell_density(end)/1e6;

    %% results
    fprintf('\n=== BIOFILM ATTACHMENT RESULTS ===\n');
    fprintf('Catheter: %.1f cm long, %.1f mm diameter\n', catheter_length, catheter_diameter*10);
    fprintf('Conditions: %.1f°C, %.1f mL/min flow\n', T, urine_flow_rate);
    
    fprintf('\nFormation Rates:\n');
    fprintf('Cell attachment rate: %.3f million cells/cm²/s\n', formation_rate_cells);
    fprintf('                    = %.1f thousand cells/cm²/s\n', formation_rate_cells * 1000);
    
    if formation_rate_time < inf
        fprintf('Time to 30%% coverage: %.0f seconds (%.1f minutes)\n', ...
                formation_rate_time, formation_rate_time/60);
    end
    
    fprintf('\nFinal State:\n');
    fprintf('Average cell density: %.1f million cells/cm²\n', avg_cell_density);
    fprintf('Maximum cell density: %.1f million cells/cm²\n', max(cell_density(:,end))/1e6);
    fprintf('Average coverage: %.1f%%\n', avg_biofilm(end)*100);
end