function [CamPos_ST] = get_Camera_Position_Struct(video_filename)
% Retrieve camera position struct.
% Use mediainfo application to retrieve video metadata information.
%   get_mediainfo_information was designed expecifically for DJI 
%   Phamtom 3 videos, and should be compatible with any DJI metadata 
%   video information.
%   
%   Note: MediaInfo is not part of CopterCurrents, and must be installed manually.
%         Media info can be found in the following links:
%
%   Media info used: MediaInfoLib - v0.7.82
%   https://mediaarea.net/en/MediaInfo/Download/Ubuntu
%   https://mediaarea.net/en/MediaInfo/Download/Windows
%
%   In the case that MediaInfo can not be installed, a work around solution 
%   will be edit the lines 108-115 manually for every video (Not recomended).  
%
%   Input:
%     video_filename: video file name
%
%   Output:
%     LON: metadata video longitude in degrees
%     LAT: metadata video latitude in degrees
%     Height: metadata  video altitude in degrees
%     timestamp: metadata Encoded date timestamp
%     mediainfo_string: raw metadata string
%     yaw: metadata yaw in degrees (heading)
%     pitch: metadata pitch in degrees 
%     roll: metadata roll in degrees 
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

command_media_info = ['mediainfo ' video_filename];

% read mediainfo information 

[status,mediainfo_string] = system(command_media_info);

if status == 0
    % get lon lat

    [C,matches] = strsplit(mediainfo_string,{'\n'});
    % Index = find(contains(C,'©xyz'));
    Index = find(contains(C,'xyz'));
    [C2,matches] = strsplit(C{Index},{':','+','-'});

    % translate coordenantes to numbers
    LON = str2double([matches{3} C2{4}]);
    LAT = str2double([matches{2} C2{3}]);
    Height = str2double([matches{4} C2{5}]);

    % search for 'Encoded date'
    Index = find(contains(C,'Encoded date'));
    [C3,matches] = strsplit(C{Index(1)},{': UTC '});
    timestamp = datenum(C3{2});

    % get orientation
    % [C,matches] = strsplit(mediainfo_string,{'\n'});
    % Index = find(contains(C,'©gyw '));
    Index = find(contains(C,'gyw '));
    [C4,matches] = strsplit(C{Index(1)},{':'});
    yaw = str2double(C4{2});

    % get pitch
    % Index = find(contains(C,'©gpt '));
    Index = find(contains(C,'gpt '));
    [C5,matches] = strsplit(C{Index(1)},{':'});
    pitch = str2double(C5{2});
    
    % get roll
    % Index = find(contains(C,'©grl '));
    Index = find(contains(C,'grl '));
    [C6,matches] = strsplit(C{Index(1)},{':'});
    roll = str2double(C6{2});
    
    
else
    % probably media info is not installed
    %  add values manually
    disp('Install Mediainfo: https://mediaarea.net/en/MediaInfo/Download/Ubuntu');
    disp('Install Mediainfo: https://mediaarea.net/en/MediaInfo/Download/Windows');
    warning('get_mediainfo_information: Mediainfo does not return valid data');

    LON = [];
    LAT = [];
    Height = [];
    timestamp = [];
    mediainfo_string = [];
    yaw = [];
    pitch = [];
    roll = [];

end

% save data in Camera Position Struct
CamPos_ST = struct('LONGITUDE',LON,'LATITUDE',LAT,'Height',Height,...
                   'timestamp',timestamp,'yaw',yaw,'pitch',pitch,...
                   'roll',roll,'extra',mediainfo_string);


end

