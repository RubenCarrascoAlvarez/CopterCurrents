function [dispersion] = get_dispersion_on_2D_fit(SNR_2D,Ux_2D,Uy_2D,lin_ind,dispersion_conf_area_mps)

% normalize SNR_2D

Ux = Ux_2D(lin_ind);
Uy = Uy_2D(lin_ind);

diff_Ux = Ux_2D - Ux;
diff_Uy = Uy_2D - Uy;

distance_Uabs = sqrt(diff_Ux.^2 + diff_Uy.^2);

% data inside confidence area
mask_conf = distance_Uabs<= dispersion_conf_area_mps;

SNR_inside_conf_area = SNR_2D(mask_conf);
SNR_outside_conf_area = SNR_2D(~mask_conf);

dispersion = sum(SNR_inside_conf_area(:)) / sum(SNR_outside_conf_area(:));

% power the SNR_2D U distance
% SNR_2D_powdist = distance_Uabs2 .* SNR_2D;
% dispersion = SNR_2D(lin_ind) /nansum(SNR_2D_powdist(:));

end

