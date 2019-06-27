function [w] = retrieve_w_from_dispersion_relation(Kx,Ky,Ux,Uy,water_depth,Domain)
% retrieve w from fundamental mode dispersion relation.
%
% w =  sqrt(9.80665 .* K .* tanh(K .* water_depth)) + (Kx .* Ux) + (Ky .* Uy);
% 
% Based on DispShell.m from Jose C. Nieto-Borge (28/06/2011)  
%      
%
% Input
%   Kx: wave number (vector or scalar) [rad/meter];
%   Ky: wave number (vector or scalar) [rad/meter];
%   Ux: Ux [m/s]. If Kx and Ky are vectors, then Ux must be a scalar;
%   Uy: Uy [m/s]. If Kx and Ky are vectors, then Uy must be a scalar;
%   water_depth: water depth in meters. If Kx and Ky are vectors, 
%                then water_depth must be a scalar;
%   Domain: +1 (positive domain) or -1 negative domain.
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

K = sqrt(Kx.^2 +  Ky.^2);
w =  (Domain .* sqrt(9.80665 .* K .* tanh(K .* water_depth))) + (Kx .* Ux) + (Ky .* Uy);

end

