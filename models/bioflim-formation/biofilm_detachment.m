
function [L] = biofilm_detachment(biofilm_density, t, L, dt, growth_velocity)
    dLdt = biofilm_detachment_rate(t, L, biofilm_density, growth_velocity);
    L = L + dLdt * dt;
    L(L<0) = 0;
    L(isnan(L)) = 0;
end

function [dLdt] = biofilm_detachment_rate(time, L, biofilm_density, growth_velocity)
    detachment_coefficient = 100 / 24 / ((100)^(-0.105)); % gCOD^(-0.965) cm^(-0.105)) h^(-1)
    shear_stress = 129600; % kg/cm h^2 
    EPS = 0.8;

    surface_detachment_rate = (detachment_coefficient * shear_stress .* (L.^2)) ./ ((biofilm_density .* EPS).^0.035);

    SOLR = 1/24000; % gCOD.cm^-2. h^-1
    k_p = 1e-6;
    detachment_period = k_p/SOLR;
    constant_volume_detachment_rate = 1/3; % per h, 8 per day

    time_dependent_detachment_coefficient = constant_volume_detachment_rate * cos((pi*time/detachment_period)^(4000 * detachment_period^2)); 
    basal_layer_thickness = 1e-4;
    volume_detachment_rate = time_dependent_detachment_coefficient .* (L - basal_layer_thickness); 

    dLdt = growth_velocity - surface_detachment_rate - volume_detachment_rate; 
    dLdt = dLdt(:);
end