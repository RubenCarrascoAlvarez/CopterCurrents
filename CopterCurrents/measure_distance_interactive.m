function [distance, distance_x,distance_y, initial_point,end_point] = measure_distance_interactive(X_2D,Y_2D,IMG_2D,maximum_n_measures)
% Tool to measure distances in a interactive way. 
%
%   Follow the istructions on the figure.
%   1) Choose area by zoom and press spacer.
%   2) Rigth click mouse in the intial point.
%   3) Rigth click mouse in the end point.
%   4) distance is displayed:
%      option 1: close figure + press Spacer => to exit 
%      option 2: press Spacer => to try again until  maximum_n_measures is
%      reaceched.
%
%   Inputs:
%     X_2D: 2 dimensional X grid 
%     Y_2D: 2 dimensional Y grid 
%     IMG_2D: grey scale image corresponding to X_2D and Y_2D.
%     maximum_n_measures: maximum tries [inf]
%
%   Outputs:
%     distance: distance retrieved sqrt(x^2+y^2), only for last try.
%     distance_x: distance in X direction (Horizontal) retrieved, 
%                 only for last try.
%     distance_y: distance in Y direction (Vertical) retrieved, 
%                 only for last try.
%     initial_point: initial mouse click coordenates,only for last try.
%     end_point: end mouse click coordenates,only for last try.
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


if exist('maximum_n_measures','var') == 0
    maximum_n_measures = inf;
end

distance = NaN;
distance_x = NaN;
distance_y = NaN;
% measure_again_flag = 1;

h_fig = figure('units','normalized','outerposition',[0 0 1 1]);
% h_fig = figure;
h_im = pcolor(X_2D,Y_2D,IMG_2D); shading flat;
xlabel('x distance');
ylabel('y distance');
h_cb = colormap(gray);
axis xy equal tight;

while ishghandle(h_fig) && maximum_n_measures> 0
    
    title('Adjust zoom, and press the SPACER continue.');
    pause;
    
    title('Left click mouse to choose initial point.');
    [x_initial,y_initial] = ginput(1);
    initial_point = [x_initial y_initial];
    
    title('Left click mouse to choose end point.');
    [x_end,y_end] = ginput(1);
    end_point = [x_end y_end];
    
    distance_x =  x_end - x_initial;
    distance_y =  y_end - y_initial;
    distance = sqrt(distance_x.^2  + distance_y.^2);
    
    if maximum_n_measures ~= 1
        title({['distance: ' num2str(distance,'%.3f')   '   distanceX: ' ...
                num2str(distance_x,'%.3f')  '   distanceY: ' ...
                num2str(distance_y,'%.3f')],...
                ['Press SPACER to continue or  CLOSE FIGURE + SPACER to finish.']});
    else
        title({['distance: ' num2str(distance,'%.3f')   '   distanceX: ' ...
                num2str(distance_x,'%.3f')  '   distanceY: ' ...
                num2str(distance_y,'%.3f')],...
                ['Press SPACER to continue']});
    end
    hold on;
    h_line = plot([x_initial x_end],[y_initial y_end],'r','linewidth',3);
    h_point = plot([x_initial x_end],[y_initial y_end],'+b');
    pause;
    
    % disable plots
    try
        set(h_line,'Visible','off');
        set(h_point,'Visible','off');
    end
    
    maximum_n_measures = maximum_n_measures - 1;
end

% close figure
if ishandle(h_fig)
    close(h_fig)
end

end

