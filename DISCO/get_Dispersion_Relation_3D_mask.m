function [DS_3D_mask] = get_Dispersion_Relation_3D_mask(Spectrum, water_depth, Ux, Uy, w_width)
% Get dispersion realtion mask from according to water_depth, Ux, Uy and w_width
% Based in Filter3DMasking.m from Jose C. Nieto-Borge (28/06/2011)  
%
% Input:
%   Spectrum: Spectrum structure from retrieve_power_spectrum
%   water_depth: water depth in meters
%   Ux: current in x direction (dim 1) [m/s];
%   Uy: current in y direction (dim 2) [m/s];
%   w_width: filter width in [rad/sec]
%
% Output: 
%   DS_3D_mask: logic mask same size tham Spectrum.power_Spectrum, where is
%   set to 1 the areas corresponding to the dispersion relation
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


% create output mask matrix
DS_3D_mask = false(size(Spectrum.W_3D));

% apply dispersion realation in K 2D space
base_DispShel_w = retrieve_w_from_dispersion_relation(Spectrum.Kx_3D(:,:,1),Spectrum.Ky_3D(:,:,1),Ux,Uy,water_depth,1);

% retrieve K 2D bouundaries according filter width
wMin_2D = base_DispShel_w  - w_width;
wMax_2D = base_DispShel_w  + w_width;

% create 3D matrix w boundaries
wMin_3D = repmat(wMin_2D,[1 1 size(DS_3D_mask,3)]);
wMax_3D = repmat(wMax_2D,[1 1 size(DS_3D_mask,3)]);


% mask values to 1
DS_3D_mask((Spectrum.W_3D >= wMin_3D) & (Spectrum.W_3D <= wMax_3D)) = 1;

% process negative part of the spectra 
if any(base_DispShel_w(:)<0) 
    % retrieve negative domain
    base_DispShel_w = retrieve_w_from_dispersion_relation(Spectrum.Kx_3D(:,:,1),Spectrum.Ky_3D(:,:,1),Ux,Uy,water_depth,-1);

    % retrieve K 2D bouundaries according filter width
    wMin_2D = base_DispShel_w  - w_width;
    wMax_2D = base_DispShel_w  + w_width;

    % create 3D matrix w boundaries
    wMin_3D = repmat(wMin_2D,[1 1 size(DS_3D_mask,3)]);
    wMax_3D = repmat(wMax_2D,[1 1 size(DS_3D_mask,3)]);


    % mask values to 1
    DS_3D_mask((Spectrum.W_3D >= wMin_3D) & (Spectrum.W_3D <= wMax_3D)) = 1;

%     cw = ceil(size(Spectrum.Kx_3D,3)/2);    
%     figure;pcolor(Spectrum.Kx_3D(:,:,1),Spectrum.Ky_3D(:,:,1),double(DS_3D_mask(:,:,cw))); shading flat;
%     cy = ceil(size(Spectrum.Kx_3D,2)/2);
%     figure;imagesc(double(squeeze(DS_3D_mask(:,cy,:)))); 
end


end

