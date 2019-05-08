function [img_monotonic_grid,conv_monotonic_grid_ST,grid_X_RW_2d,grid_Y_RW_2d] ...
    = georeferenceDJIFrame_byCaltech_EQ_v2(img,altitude,pitch,roll,DISCO_CamCalib)
% Get X and Y Real world coordeantes grid from a given Image using 
% DISCO_CamCalib note X and Y grid monotonic grid, and the reference is located in
% the center of the image (DISCO_CamCalib.cc, principal point coordenates)
%
%   Note: this function uses the 'Camera Calibration Toolbox for Matlab'  
%   an external libary, written by Jean-Yves Bouguet,
%   jean-yves.bouguet@intel.com. This library is not part of DISCO, and it
%   is available in:
%
%   http://www.vision.caltech.edu/bouguetj/calib_doc/
%   
%   The 'Camera Calibration Toolbox for Matlab' is compatible with OpenCV
%   camera calibration. This is the main reason to include this function in
%   DISCO.
%   
%   If you want to avoid external DISCO libraries, then use the function
%   'gridDJIFrame_nadir.m' instead.
%
% 
% Intput:
%   img: grey scale image.
%   altitude: distance from the camera to the water surface in meters.
%   pitch: camera pitch in degrees. (-90 degrees for nadir adquisition)
%   roll: camera roll in degrees. (0 degrees for nadir adquisition)
%   DISCO_CamCalib: calibration structure from caltech toolbox or openCV
%
% Output:
%    img_monotonic_grid: image converted to a monotonic grid. The first
%       img_monotonic_grid(x,y)
%    conv_monotonic_grid_ST: conversion to monotonic grid structure
%         conv_monotonic_grid_ST.LinearInd_grid: linear index used for 
%                                                conversion. 
%         conv_monotonic_grid_ST.LinearInd_img: corresponding linear index
%                    to 'LinearInd_grid' in the original grey scale image.
%         conv_monotonic_grid_ST.grid_d1_RW_2d: 1 dimension Real Wold 2d grid. 
%         conv_monotonic_grid_ST.grid_d2_RW_2d: 2 dimension Real Wold 2d grid.  
%         conv_monotonic_grid_ST.size_original_image: size of the original image.
%         conv_monotonic_grid_ST.dxdy: outout grid resoultion.
%    grid_X_RW_2d: X dimension Real Wold 2d grid (x,y). 
%    grid_Y_RW_2d: Y dimension Real Wold 2d grid (x,y). 
%
%   note: use 'apply_conv_monotonic_grid_ST' to retrieve the image in (x,y)
%         format, with axis_flag = 0 
%   [X_eq,Y_eq,IMG_eq] = apply_conv_monotonic_grid_ST(...
%                      conv_monotonic_grid_ST, img_gray,axis_flag);
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


% check image size
if size(img,1) ~= DISCO_CamCalib.ny
    error('dim1: Y axis (vertical) does not match with DISCO_CamCalib.ny');
end

if size(img,2) ~= DISCO_CamCalib.nx
    error('dim2: X axis (horizontal) does not match with DISCO_CamCalib.nx');
end

if ndims(img)>3
    error('The img must be only 2 dimesions (grey sacale image) ');
end


% create camera matrix
% KK = [DISCO_CamCalib.fc(1)  DISCO_CamCalib.fc(1)*DISCO_CamCalib.alpha_c  DISCO_CamCalib.cc(1);...
%            0                DISCO_CamCalib.fc(2)                         DISCO_CamCalib.cc(2);...
%            0                            0                                           1];

% Rotation vector
% https://developer.dji.com/mobile-sdk/documentation/introduction/flightController_concepts.html
% Roll => X
% Pitch => Y
% Yaw => Z
om = deg2rad([roll pitch 0]);
% om = deg2rad([0 0 pitch]);
if roll ~= 0
   warning('georeferenceDJIFrame_byCaltech_v2: expect roll = 0.') 
end

% get Translation vector
% T = [0; 0 ; altitude];
T = [0; 0 ; 0];

% note: X and Y are exchange to match with IMG sequence notation. 
% (X => horizontal)
% (Y => vertical)
% point coordenates [X Y Z] => project_points3

% get X limits 
RatAlt = 1;
X_lim = [NaN NaN];
while(any(isnan(X_lim)))
    vec = -RatAlt*altitude: altitude/1000:RatAlt*altitude; % create vector [X 0 altitude];
    X = cat(2,vec(:),zeros(numel(vec),1),repmat(altitude,[numel(vec) 1]))';
    [xp] = project_points3(X,om,T,DISCO_CamCalib.fc,DISCO_CamCalib.cc,DISCO_CamCalib.kc,DISCO_CamCalib.alpha_c);
    % value of Y where the pixel coordenate is 1
    % note: X is conected with the dim1 of xp
    X_lim = interp1(xp(1,:),vec,[1 DISCO_CamCalib.nx]-1,'linear');
    RatAlt = RatAlt + 1;
end
X_lim = sort(X_lim);


% point coordenates [X Y Z] => project_points3
% get Y limits 
RatAlt = 1;
Y_lim = [NaN NaN];
while(any(isnan(Y_lim)))
    vec = -RatAlt*altitude: altitude/1000 :RatAlt*altitude; % crate vector [0 Y altitude];
    Y = cat(2,zeros(numel(vec),1),vec(:),repmat(altitude,[numel(vec) 1]))';
    [xp] = project_points3(Y,om,T,DISCO_CamCalib.fc,DISCO_CamCalib.cc,DISCO_CamCalib.kc,DISCO_CamCalib.alpha_c);
    % value of Y where the pixel coordenate is 1
    % note: Y is conected with the dim2 of xp
    Y_lim = interp1(xp(2,:),vec,[1 DISCO_CamCalib.ny]-1,'linear');
end
Y_lim = sort(Y_lim);

% create monotonic and ascending grid (Real World Coordenates)
offsetRatio = 1.0;
grid_Y_RW_1d = linspace(Y_lim(1)*offsetRatio ,Y_lim(2)*offsetRatio ,DISCO_CamCalib.ny);
dxdy = grid_Y_RW_1d(2) - grid_Y_RW_1d(1);
grid_X_RW_1d = X_lim(1)*offsetRatio : dxdy : X_lim(2)*offsetRatio;


% create equidistant and monotonic grid matrix
% [grid_X_RW_2d ,grid_Y_RW_2d ] = meshgrid(grid_X_RW_1d ,grid_Y_RW_1d);
[grid_X_RW_2d ,grid_Y_RW_2d ] = ndgrid(grid_X_RW_1d ,grid_Y_RW_1d);


% create inputs
X = cat(2,grid_X_RW_2d(:),grid_Y_RW_2d(:),ones(numel(grid_X_RW_2d),1)*altitude)';

% use caltech funcion (v2) project_points2 x y coordenates to pixels
%  [xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk] = project_points2(X,om,T,f,c,k,alpha)
[xp] = project_points2(X,om,T,...
    DISCO_CamCalib.fc,DISCO_CamCalib.cc,DISCO_CamCalib.kc,DISCO_CamCalib.alpha_c);

% [xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk] = project_points(X,om,T,f,c,k)
% [xp] = project_points(X,om,T,DISCO_CamCalib.fc,DISCO_CamCalib.cc,DISCO_CamCalib.kc);

% reshape (X,Y)
pixels_2D_X =  reshape(xp(1,:),size(grid_X_RW_2d)); % 
pixels_2D_Y  = reshape(xp(2,:),size(grid_X_RW_2d));

% get nearest pixels
pixels_2D_X = round(pixels_2D_X)+1;
pixels_2D_Y = round(pixels_2D_Y)+1;

%  set to nan boundaries out of 'img'
pixels_2D_X(pixels_2D_X<1 | pixels_2D_X >DISCO_CamCalib.nx) = NaN;
pixels_2D_Y(pixels_2D_Y<1 | pixels_2D_Y >DISCO_CamCalib.ny) = NaN;

% flip Y axis (due to the image is fliped by video reader)
pixels_2D_Y = flipud(pixels_2D_Y);

% get linear index corresponding to the orinal image
% LinearInd_img = sub2ind(size(img),pixels_2D_Y(:),pixels_2D_X(:));
LinearInd_img = sub2ind(size(img),pixels_2D_Y(:),pixels_2D_X(:));

% get linear index corresponded to the new grids [grid_d1_RW_2d,grid_d2_RW_2d]
[grid_X_RW_2d_pix,grid_Y_RW_2d_pix] = meshgrid(1:length(grid_X_RW_1d),1:length(grid_Y_RW_1d));
LinearInd_grid = sub2ind(size(grid_X_RW_2d_pix),grid_Y_RW_2d_pix(:),grid_X_RW_2d_pix(:));


% keep only valid index 
ind2use = isfinite(LinearInd_img);
LinearInd_grid = LinearInd_grid(ind2use);
LinearInd_img = LinearInd_img(ind2use);

% apply conversion
img_monotonic_grid = nan(size(pixels_2D_X));
img_monotonic_grid(LinearInd_grid) = img(LinearInd_img);

% store conversion structure
size_original_image = size(img);
grid_d1_RW_2d = grid_X_RW_2d;
grid_d2_RW_2d = grid_Y_RW_2d;

conv_monotonic_grid_ST =  struct('LinearInd_grid',LinearInd_grid,...
    'LinearInd_img',LinearInd_img,'grid_d1_RW_2d',grid_d1_RW_2d,...
    'grid_d2_RW_2d',grid_d2_RW_2d,'size_original_image',...
    size_original_image,'dxdy',dxdy);


% % plot
% figure;
% pcolor(grid_X_RW_2d ,grid_Y_RW_2d,img_monotonic_grid);
% shading flat;
% axis xy equal tight;
% colormap(gray);


end