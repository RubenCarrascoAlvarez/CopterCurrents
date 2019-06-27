function [Spectrum] = retrieve_power_spectrum(IMG_3D,dx,dy,dt,K_limits,W_limits)
% Retrieve 3D structure with power spectrum. The power spectrum is cutted
% according to K_limits, W_limits and normalized where
%   sum(Spectrum.power_Spectrum(:)) = 1 
% 
%
%   Input
%    IMG_3D: 3D array with image sequence (X,Y,Time)
%    dx: resolution in x (dim 1) [meters] 
%    dy: resolution in y (dim 2) [meters] 
%    dt: resolution in time (dim 3) [seconds] 
%    K_limits: [min max] Wave-number interval to retrieve [rad/m]
%    W_limits: [min max] frequency interval to retrieve [rad/sec]
%
%   Output
%    Spectrum.power_Spectrum: (Kx,Ky,W) power Spectrum
%    Spectrum.Kx_3D: 3D Kx grid corresponding to Spectrum.power_Spectrum [rad/m]
%    Spectrum.Ky_3D: 3D Ky grid corresponding to Spectrum.power_Spectrum [rad/m]
%    Spectrum.W_3D: 3D W grid corresponding to Spectrum.power_Spectrum [rad/sec]
%    Spectrum.dKx: Kx resolution [rad/m]
%    Spectrum.dKy: Ky resolution [rad/m]
%    Spectrum.dW: W resolution [rad/sec]
%    Spectrum.Kx_orig_limits : orininal Kx vector limits [rad/m]
%    Spectrum.Ky_orig_limits : orininal Ky vector limits [rad/m]
%    Spectrum.W_orig_limits: orininal w vector limits [rad/sec] 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of CopterCurrents
%
% CopterCurrents has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% CopterCurrents is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CopterCurrents is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CopterCurrents.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if exist('K_limits','var') ==0 || isempty(K_limits)
    K_limits = [-inf inf];
end

if exist('W_limits','var') ==0 || isempty(W_limits)
    W_limits = [-inf inf];
end

% check limits
if any(K_limits<0) 
    warning('retrieve_positive_spectrum: K_limits smaller than 0');
end

if any(W_limits<0) 
    warning('retrieve_positive_spectrum: W_limits smaller than 0');
end

% get Kx Ky abd w axis
[Nx, Ny, Nt] = size(IMG_3D);
Kx = (2 * pi * 1/dx/Nx) * [-ceil((Nx-1)/2): floor((Nx-1)/2)];
Ky = (2 * pi * 1/dy/Ny) * [-ceil((Ny-1)/2): floor((Ny-1)/2)];
w  = (2 * pi * 1/dt/Nt) * [-ceil((Nt-1)/2): floor((Nt-1)/2)];

dKx = Kx(2)-Kx(1);
dKy = Ky(2)-Ky(1);
dW = w(2)-w(1);
Norm = dKx * dKy * dW;

% get power spectra
Spectrum_raw = fftshift(fftn(IMG_3D)/numel(IMG_3D));
power_Spectrum = abs(Spectrum_raw/Norm).^2; % normalization not neccesary

% get index boundaries
 % K domain
ind_x = abs(Kx) <= K_limits(2);
ind_y = abs(Ky) <= K_limits(2);
 % W domain
ind_w = w  >= W_limits(1) &  w <= W_limits(2);


% create 3D structure
% [Kx_2D,Ky_2D] = ndgrid(Kx(ind_x),Ky(ind_y));
[Kx_3D,Ky_3D,W_3D] = ndgrid(Kx(ind_x),Ky(ind_y),w(ind_w));
power_Spectrum_cut = power_Spectrum(ind_x,ind_y,ind_w);

% set to nan values out of K_limits
K_3D = sqrt(Kx_3D.^2 + Ky_3D.^2);
power_Spectrum_cut(K_3D<K_limits(1) | K_3D>K_limits(2)) = NaN; 

% get original spectrum limits
Kx_orig_limits = [Kx(1) Kx(end)];
Ky_orig_limits = [Ky(1) Ky(end)];
W_orig_limits = [w(1) w(end)];

% normalize spectra (sum(power_Spectrum_cut(:)) = 1)
power_Spectrum_cut = power_Spectrum_cut / nansum(power_Spectrum_cut(:));

% create output structure

Spectrum = struct('power_Spectrum',power_Spectrum_cut,...
           'Kx_3D',Kx_3D,'Ky_3D',Ky_3D,'W_3D',W_3D,...
           'dKx',dKx,'dKy',dKy,'dW',dW,...
           'Kx_orig_limits',Kx_orig_limits,'Ky_orig_limits',...
            Ky_orig_limits,'W_orig_limits',W_orig_limits);

% ind_w_pos = w>=0;
% % ind_w_pos = w>=W_limits(1);
% W_full_positive = w(ind_w_pos);
% power_Spectrum_full_positive = power_Spectrum(:,:,ind_w_pos);
% 
% % % Normalize spectra to sum = 1
% power_Spectrum_cut = power_Spectrum_cut / nansum(power_Spectrum_cut(:));
% power_Spectrum_full_positive = power_Spectrum_full_positive / nansum(power_Spectrum_full_positive(:));
% 
% Spectrum = struct('power_Spectrum',power_Spectrum_cut,...
%            'Kx_3D',Kx_3D,'Ky_3D',Ky_3D,'W_3D',W_3D,...
%            'dKx',dKx,'dKy',dKy,'dW',dW,...
%            'Kx_full',Kx,'Ky_full',...
%             Ky,'W_full_positive',W_full_positive,...
%            'power_Spectrum_full_positive',power_Spectrum_full_positive);


end

