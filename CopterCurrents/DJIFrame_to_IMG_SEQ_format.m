function [IMG_SEQ_frame] = DJIFrame_to_IMG_SEQ_format(DJI_frame)
%   Transform DJI frame to IMG_SEQ_frame format. Two transformation are
%   needed:
%   - flip Y dimension
%   - permute dimension 1 by 2.  
%   Input:
%     DJI_frame: grey scale raw 2D Intensity image (Ny,Nx), extrated from   
%          DII video:
%          dim1 => Y (vertical) data. Monotonic descending.
%          dim2 => X (Horizontal) data. Monotonic ascending.
%
%
%   Output:
%     IMG_SEQ_frame: 2D matrix (Nx,Ny).
%          dim1 => X (Horizontal) data. Monotonic ascending.
%          dim2 => Y (vertical) data.  Monotonic ascending.
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

% perform flip Y axis, and permute 1 and 2 dimension.
IMG_SEQ_frame = permute(flipud(DJI_frame),[2 1]);

end

