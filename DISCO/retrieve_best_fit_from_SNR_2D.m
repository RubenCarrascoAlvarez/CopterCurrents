function [Ux,Uy,SNR_max,lin_ind] = retrieve_best_fit_from_SNR_2D(SNR_2D,Ux_2D,Uy_2D,method,display)
% Retrieve best fit from SNR 2D maps.
%
%   Input:
%       SNR_2D: signal to noise ration matrix 
%       Ux_2D: 2D Ux matrix corresponding to SNR_2D
%       Uy_2D: 2D Uy matrix corresponding to SNR_2D
%       method: method to find the best fit
%               0 => maximum Signal to Noise Ratio
%               1 => Not implemented yet.
%
%   Output:
%       Ux: Best fit for Ux 
%       Uy: Best fit for Uy 
%       SNR_max: SNR corresponding to the best fit.
%       lin_ind: 2D matrix linear index corresponding to the best fit.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of DISCO
%
% DISCO has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% DISCO is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% DISCO is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DISCO.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method == 0 % find maximum SNR
    
    % find maximum SNR
    [SNR_max,lin_ind] = max(SNR_2D(:));
    
    % retrieve cor
    Ux = Ux_2D(lin_ind);
    Uy = Uy_2D(lin_ind);
    
elseif method == 1 % add more methods (TO DO)
    Ux = [];
    Uy = [];
    
    error(['retrieve_best_fit_from_SNR_2D: unknown method ' num2str(method)]);
end

% display = 1;
% display = 0;

if display == 1
    
   offset_x =  (Ux_2D(2,1) - Ux_2D(1,1))/2;
   offset_y =  (Uy_2D(1,2) - Uy_2D(1,1))/2;
   
   figure;
   pcolor( Ux_2D - offset_x, Uy_2D - offset_y,SNR_2D ); 
   shading flat;
   axis xy equal tight;
   hc = colorbar;
   ylabel('SNR');
   hold on;
   plot(Ux,Uy,'*k');
   xlabel('Ux [m/s]');
   ylabel('Uy [m/s]');
   title({['Ux: ' num2str(Ux) ' m/s'],['Uy: ' num2str(Uy) ' m/s']})
    
end

end

