function [GeoMode] = getGeoMode_from_CopterCurrents_CamCalib(CopterCurrents_CamCalib)
% Get the type of georeference to be used in the CopterCurrents_CamCalib structure
%
%   Input:
%       CopterCurrents_CamCalib: CopterCurrents Camera calibration structure.
%
%   Output:
%       GeoMode: 0 => Caltech or OpenCV georeferencing method.
%       GeoMode: 1 => FOV georeferencing method.

% GeoMode = 0 => 
%           CopterCurrents_CamCalib.fc: focal length in pixels, 2x1 vector.
%           CopterCurrents_CamCalib.cc: principal point coordinates, 2x1 vector.
%           CopterCurrents_CamCalib.Skew coefficient: The skew coefficient defining
%               the angle  between the x and y pixel axes, 1x1 scalar, 
%               normally 0.
%       CopterCurrents_CamCalib.Distortions: Image distortion coefficients
%           (radial and tangential distortions), 5x1 vector where:
%               kc(1) = k1 (radial distortion)
%               kc(2) = k2 (radial distortion)
%               kc(3) = p1 (tangencial distortion)
%               kc(4) = p2 (tangencial distortion)
%               kc(5) = k3 (radial distortion)
%       CopterCurrents_CamCalib.nx: number of sensor pixels in horizontal direction. 
%       CopterCurrents_CamCalib.ny: number of sensor pixels in vertical direction.
%       CopterCurrents_CamCalib.camera_offset_Z: distance in meters between the
%               camera and drone GPS antenna.
%       CopterCurrents_CamCalib.source: Calibration source used.
%
% GeoMode = 1 => 
%    CopterCurrents_CamCalib: camera configuration structure
%         CopterCurrents_CamCalib.fov_diag: diagonal field of view angle
%         CopterCurrents_CamCalib.fov_x : x (horizontal) field of view angle;
%         CopterCurrents_CamCalib.fov_y : y (vertical) field of view angle;
%         CopterCurrents_CamCalib.size_X : frame pixel size in x (horizontal);
%         CopterCurrents_CamCalib.size_Y : frame pixel size in y (vertical);
%         CopterCurrents_CamCalib.camera_offset_Z = vertical (Z) distance in meter
%                                          between GPS and camera;
%         CopterCurrents_CamCalib.source = 'Manual_FOV_calibration';
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




% check fov_diag field to know which georeferencing method should be used.
if isfield(CopterCurrents_CamCalib,'fov_diag')
    GeoMode = 1;
else
    GeoMode = 0;
    disp(CopterCurrents_CamCalib.source);
end

disp(['Georeferencing Mode: ' num2str(GeoMode)]);

end

