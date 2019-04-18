% add DISCO folders and subfolders to matlab path
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


%%%%%%%%%%%%%%%%%%%%%%%%%
% Add to the matlab path
%%%%%%%%%%%%%%%%%%%%%%%%%

% get all the folders in the DISCO path
DISCO_main_path = [pwd filesep];
paths2add_str = genpath(DISCO_main_path);

% list of folder to add to the path
% valid_list_folders = {'external_libraries'};
valid_list_folders = {'external_libraries','test_scripts'};

% check if extra folder is need
if exist('nanmean','file') ~= 2 || exist('nansum','file') ~= 2
    valid_list_folders = cat(2,valid_list_folders,'extra');
end
 
% keep only paths included in valid_list_folders
str_delimited = split(paths2add_str,':');
n_match_found = zeros(size(str_delimited));

for i1 = 1:length(valid_list_folders) 
    fstr = strfind(str_delimited,valid_list_folders{i1});
    n_match_found = n_match_found + double(~cellfun(@isempty,fstr));
end

% remove folders without matches
str_delimited(n_match_found == 0) = [];

% generate new genpath string
paths2add_str_clean = [DISCO_main_path ':'];
for i1 = 1:length(str_delimited) 
    paths2add_str_clean = [paths2add_str_clean str_delimited{i1} ':'];
end

addpath(paths2add_str_clean);

%%%%%%%%%%%%%%%%%
% check options
%%%%%%%%%%%%%%%%%
% clear displayed messages
clc

% mandatory libraries

% show matlab version and toolboxes
ver

disp(' ');

% check mediainfo
% [status,result] = system('mediainfo');
[status,result] = system('mediainfo --Version');
disp(result);

if status ~= 255
   warning('Mediainfo probably not installed.') 
end

% check nanmean function availability
if exist('nanmean','file') == 2 
    disp('nanmean function found.');
else
    warning('nanmean function NOT found.');
end
    
    
% check nasum function availability
if exist('nansum','file') == 2 
    disp('nansum  function found.');
else
    warning('nansum  function NOT found.');
end
disp(' ');


% optional libraries
disp('----------------------------------------------');
disp('Checking optional libraries:');
disp('----------------------------------------------');

% deg2utm
if exist('deg2utm','file') == 2
    disp('deg2utm: AVAILABLE (Optional)'); 
else
    disp('deg2utm: NOT FOUND (Optional)'); 
end

% Camera Calibration Toolbox
if exist('calib_gui','file') == 2
    disp('Camera Calibration Toolbox: AVAILABLE (Optional)'); 
else
    disp('Camera Calibration Toolbox: NOT FOUND (Optional)'); 
end
