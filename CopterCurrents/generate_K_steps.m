function [K_steps] = generate_K_steps(wavelength_min,wavelength_max,K_width,K_distance)

% wavelength_min: minimum wavelenght in meters
% wavelength_max: maximum wavelenght in meters
% K_width: bandwidth of the K_steps in wavenumbers [defout value 1]
% K_distance: distance between 2 K_steps in wave numbers [K_width/2]


if exist('K_width','var') == 0 || isempty(K_width)
    K_width = 1;
end

if exist('K_distance','var') == 0 || isempty(K_distance)
    K_distance = 1;
end

K_ini = 2*pi/wavelength_max;
K_end = 2*pi/wavelength_min;

step_ini = K_ini:K_distance:K_end;
step_end = step_ini + K_width;
K_steps = cat(1,step_ini,step_end);

% remove steps out of wavelength_min wavelength_max boundaries
ind2del_K_steps =  any(K_steps<K_ini | K_steps>K_end,1);
K_steps(:,ind2del_K_steps) = [];

end

