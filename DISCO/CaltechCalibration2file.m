function [DISCO_CamCalib] = CaltechCalibration2file(Caltech_calibration_filename, DISCO_calibration_filename,camera_offset_Z)
%   Translate Caltech calibration files (genrated by the 'Camera 
%   Calibration Toolbox for Matlab') to a DISCO calibration camera file.
%
%   Input:
%    Caltech_calibration_filename: calibration file generated my Caltech
%    matlab calibration GUI:
%
%     (http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example.html)
%
%    DISCO_calibration_filename: DISCO camera calibration filename, where
%    the intrinsic camera parameters will be stored.
%    camera_offset_Z: distance in meters between the camera and drone 
%                     GPS antenna.
%
%   Output:
%       DISCO_CamCalib: DISCO Camera calibration structure.
%
%           DISCO_CamCalib.fc: focal length in pixels, 2x1 vector.
%           DISCO_CamCalib.cc: principal point coordinates, 2x1 vector.
%           DISCO_CamCalib.Skew coefficient: The skew coefficient defining
%               the angle  between the x and y pixel axes, 1x1 scalar, 
%               normally 0.
%       DISCO_CamCalib.Distortions: Image distortion coefficients
%           (radial and tangential distortions), 5x1 vector where:
%               kc(1) = k1 (radial distortion)
%               kc(2) = k2 (radial distortion)
%               kc(3) = p1 (tangencial distortion)
%               kc(4) = p2 (tangencial distortion)
%               kc(5) = k3 (radial distortion)
%       DISCO_CamCalib.nx: number of sensor pixels in horizontal direction. 
%       DISCO_CamCalib.ny: number of sensor pixels in vertical direction.
%       DISCO_CamCalib.camera_offset_Z: distance in meters between the
%               camera and drone GPS antenna.
%       DISCO_CamCalib.source: Calibration source used.
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

if exist('camera_offset_Z','var') == 0 || isempty(camera_offset_Z)
   camera_offset_Z = 0; 
end

% load Caltech_calibration_filename
DISCO_CamCalib = load(Caltech_calibration_filename,'fc','cc','kc','alpha_c','nx','ny');

% Add camera offset
DISCO_CamCalib.camera_offset_Z = camera_offset_Z;

% add source field
DISCO_CamCalib.source = 'Caltech_Matlab';

% save DISCO_CamCalib structure to DISCO_calibration_filename 
save(DISCO_calibration_filename,'DISCO_CamCalib');
disp(['Calibration saved: ' DISCO_calibration_filename]);

% display save structure
disp(DISCO_CamCalib);


end

