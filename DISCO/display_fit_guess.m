function display_fit_guess(IMG_SEQ,STCFIT,n_window,Ux,Uy,w_cut_dir)
%   GUI very usefull to check the parameters in STCFIT applied to the data 
%   IMG_SEQ.
%
%   Input:
%       IMG_SEQ:
%
%       STCFIT: Structure defining the current fit input parameters
%
%         STCFIT.Generic.gridX: 2D X grid (Horizontal) in meters
%         STCFIT.Generic.gridY: 2D Y grid (Vertical) in meters 
%         STCFIT.Generic.image: image example
%         STCFIT.Generic.Longitude: Longitude [deg]
%         STCFIT.Generic.Latitude:  Latitude [deg]
%         STCFIT.Generic.heading:   Angle to North [deg]
%         STCFIT.Generic.altitude:  altitude in meters
%         STCFIT.Generic.time_stamp: time stamps in datenum format
%         STCFIT.Generic.usable_data_mask_2D: usable points for fit.
%
%         STCFIT.Windows.N_fit_windows: number of windows used for the current fit.
%         STCFIT.Windows.w_corners_dim1: 3xN_fit_windows
%         STCFIT.Windows.w_corners_dim2: 3xN_fit_windows
%         STCFIT.Windows.sq_size_m: square size in meter for the current fit windows.
%         STCFIT.Windows.sq_dist_m: square distance in meter between current fit windows.
%         STCFIT.Windows.average_depth_win: average depth in meters.
% 
%         STCFIT.fit_param.Ux_FG_2D: MxN Ux first guess matrix
%         STCFIT.fit_param.Uy_FG_2D: MxN Uy first guess matrix
%         STCFIT.fit_param.Ux_SG_2D: MxN Ux second guess matrix (offset matrix)
%         STCFIT.fit_param.Uy_SG_2D: MxN Uy second guess matrix (offset matrix) 
%         STCFIT.fit_param.w_width_FG: first guess filter width in w [rad/s] 
%         STCFIT.fit_param.w_width_SG: second guess filter width in w [rad/s]
%         STCFIT.fit_param.waveLength_limits_m: [min max] wavelength to use [meters]
%         STCFIT.fit_param.wavePeriod_limits_sec: [min max] wave Period to use [seconds]
%         STCFIT.fit_param.K_limits: [min max] wave number to use [rad/m]
%         STCFIT.fit_param.W_limits: [min max] wave frequency to use [rad/sec]
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


    
if exist('Ux','var') == 0 || isempty(Ux)
    Ux = 0;
end

if exist('Uy','var') == 0 || isempty(Uy)
    Uy = 0;
end

if exist('w_cut_dir','var') == 0 || isempty(w_cut_dir)
    w_cut_dir = 90;
end

    
% cut window
IMG_SEQ_Window = double(IMG_SEQ.IMG(STCFIT.Windows.w_corners_dim1(2,n_window):STCFIT.Windows.w_corners_dim1(3,n_window),...
                                    STCFIT.Windows.w_corners_dim2(2,n_window):STCFIT.Windows.w_corners_dim2(3,n_window),:));

% get spectrum
Spectrum = retrieve_power_spectrum(IMG_SEQ_Window,IMG_SEQ.dx, IMG_SEQ.dy, IMG_SEQ.dt,STCFIT.fit_param.K_limits, STCFIT.fit_param.W_limits); 

% get parametres
water_depth = STCFIT.Windows.average_depth_win(n_window);
w_width = STCFIT.fit_param.w_width_FG;

% get K vector
K_cut = min(Spectrum.Kx_3D(:)):Spectrum.dKx:max(Spectrum.Kx_3D(:));
W_cut = squeeze(Spectrum.W_3D(1,1,:));


%%%%%%%%%%%%%%%%        
% create figure       
%%%%%%%%%%%%%%%%    
h = figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) =  subplot(2,2,1);
ax(2) =  subplot(2,2,2);
% ax(3) =  subplot(2,2,3);
ax(4) =  subplot(2,2,4);

% create text box annotations
scroll_y_size = 0.04;
scroll_y_offset = scroll_y_size * 2.0;

% ann 1
pos_ann_1 = [ax(1).Position(1) ax(4).Position(2) ax(4).Position(3)/2 scroll_y_size];
str_ann = ['filter width ' num2str(w_width) ' [rad/s]' ];
h_ann_1 = annotation('textbox', pos_ann_1, 'String', str_ann,'FontWeight','Bold',...
                     'HorizontalAlignment','center','VerticalAlignment','middle');

% ann 2
pos_ann_2 = [ax(1).Position(1) pos_ann_1(2)+scroll_y_offset ax(4).Position(3)/2 scroll_y_size];
str_ann = ['W cut direction: ' num2str(w_cut_dir) ' [deg]' ];
h_ann_2 = annotation('textbox', pos_ann_2, 'String', str_ann,'FontWeight','Bold',...
                     'HorizontalAlignment','center','VerticalAlignment','middle');

% ann 3
pos_ann_3 = [ax(1).Position(1) pos_ann_2(2)+scroll_y_offset ax(4).Position(3)/2 scroll_y_size];
str_ann = ['Uy: ' num2str(Uy) ' [m/s]' ];
h_ann_3 = annotation('textbox', pos_ann_3, 'String', str_ann,'FontWeight','Bold',...
                     'HorizontalAlignment','center','VerticalAlignment','middle');

% ann 4
pos_ann_4 = [ax(1).Position(1) pos_ann_3(2)+scroll_y_offset ax(4).Position(3)/2 scroll_y_size];
str_ann = ['Ux: ' num2str(Ux) ' [m/s]' ];
h_ann_4 = annotation('textbox', pos_ann_4, 'String', str_ann,'FontWeight','Bold',...
                     'HorizontalAlignment','center','VerticalAlignment','middle');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save plot data in userdata struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.UserData = struct('Spectrum',Spectrum,'K_cut',K_cut,'W_cut',W_cut,...
            'water_depth',water_depth,'Ux',Ux,'Uy',Uy,'w_cut_dir',w_cut_dir,...
            'w_width',w_width,'ax',ax,'h_ann_1',h_ann_1,'h_ann_2',h_ann_2,...
            'h_ann_3',h_ann_3,'h_ann_4',h_ann_4);
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    create sliding bars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create slider w_width
ypos_sc1 = [ax(1).Position(1)+(ax(4).Position(3)/2) ax(4).Position(2) ax(4).Position(3)/2 scroll_y_size];
ymin_sc1 = 0;
ymax_sc1 = 4;

clbk_sc1 = @selection_sc1;
h_sc1 = uicontrol('style','slider', ...
    'units','normalized','position',ypos_sc1, ...
    'callback',clbk_sc1,'min',ymin_sc1,'max',ymax_sc1,...
    'value',h.UserData.w_width,'SliderStep', [1/40 , 0.25]);    


% create slider w_cut_dir
ypos_sc2 = [ax(1).Position(1)+(ax(4).Position(3)/2) ypos_sc1(2)+scroll_y_offset ax(4).Position(3)/2 scroll_y_size];
ymin_sc2 = 0;
ymax_sc2 = 360;

clbk_sc2 = @selection_sc2;
h_sc2 = uicontrol('style','slider', ...
    'units','normalized','position',ypos_sc2, ...
    'callback',clbk_sc2,'min',ymin_sc2,'max',ymax_sc2,...
    'value',h.UserData.w_cut_dir,'SliderStep', [1/360 , 1/4]);    

% create slider Uy
ypos_sc3 = [ax(1).Position(1)+(ax(4).Position(3)/2) ypos_sc2(2)+scroll_y_offset ax(4).Position(3)/2 scroll_y_size];
ymin_sc3 = -3;
ymax_sc3 = 3;

clbk_sc3 = @selection_sc3;
h_sc3 = uicontrol('style','slider', ...
    'units','normalized','position',ypos_sc3, ...
    'callback',clbk_sc3,'min',ymin_sc3,'max',ymax_sc3,...
    'value',h.UserData.Uy,'SliderStep', [1/120 , 1/12]);    

% create slider Ux
ypos_sc4 = [ax(1).Position(1)+(ax(4).Position(3)/2) ypos_sc3(2)+scroll_y_offset ax(4).Position(3)/2 scroll_y_size];
ymin_sc4 = -3;
ymax_sc4 = 3;

clbk_sc4 = @selection_sc4;
h_sc4 = uicontrol('style','slider', ...
    'units','normalized','position',ypos_sc4, ...
    'callback',clbk_sc4,'min',ymin_sc4,'max',ymax_sc4,...
    'value',h.UserData.Ux,'SliderStep', [1/120 , 1/12]);    


%%%%%%%%%%%%%%%%%%
%  update figure
%%%%%%%%%%%%%%%%%%
h = update_spectra_plot(h);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Update sliding bar funtions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selection_sc2(src,event)
    event.Source.Parent.UserData.w_cut_dir = src.Value;
    disp(['New Angle: ' num2str(event.Source.Parent.UserData.w_cut_dir)]);
    update_spectra_plot(event.Source.Parent);
end

function selection_sc1(src,event)
    event.Source.Parent.UserData.w_width = src.Value;
    disp(['New Filter width: ' num2str(event.Source.Parent.UserData.w_width)]);
    update_spectra_plot(event.Source.Parent);
end

function selection_sc3(src,event)
    event.Source.Parent.UserData.Uy = src.Value;
    disp(['New Angle: ' num2str(event.Source.Parent.UserData.Uy)]);
    update_spectra_plot(event.Source.Parent);
end

function selection_sc4(src,event)
    event.Source.Parent.UserData.Ux = src.Value;
    disp(['New Angle: ' num2str(event.Source.Parent.UserData.Ux)]);
    update_spectra_plot(event.Source.Parent);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Update figure function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = update_spectra_plot(h)

    % copy userdata parameters
    Spectrum = h.UserData.Spectrum;
    K_cut = h.UserData.K_cut;
    W_cut = h.UserData.W_cut;
    Ux = h.UserData.Ux;
    Uy = h.UserData.Uy;
    water_depth = h.UserData.water_depth;
    w_cut_dir = h.UserData.w_cut_dir;
    w_width = h.UserData.w_width;
    ax = h.UserData.ax;
    
    % cut in projection direction
    Kx_cut = K_cut * sind(w_cut_dir);
    Ky_cut = K_cut * cosd(w_cut_dir);

    % get dispersion relatio mask according to Ux Uy
    [DR_3D_mask] = get_Dispersion_Relation_3D_mask(Spectrum, water_depth, Ux, Uy, w_width);

    % get power average in W
    power_Spectrum_3D_DR = Spectrum.power_Spectrum;
    power_Spectrum_3D_DR(DR_3D_mask~=1) = NaN;
    power_Spectrum_2D_DR = nansum(power_Spectrum_3D_DR,3);
    power_Spectrum_2D_DR_log = 10*log10(power_Spectrum_2D_DR);

    % display w_cut_dir
    % create interpolation structure
    F = griddedInterpolant(Spectrum.Kx_3D, Spectrum.Ky_3D, Spectrum.W_3D,Spectrum.power_Spectrum,'nearest','none');
    F_mask = griddedInterpolant(Spectrum.Kx_3D, Spectrum.Ky_3D, Spectrum.W_3D,power_Spectrum_3D_DR,'nearest','none');

    % create K,W vector to interpolate
    [Kx_cut_2D,W_cut_2D]= ndgrid(Kx_cut,W_cut);
    [Ky_cut_2D,W_cut_2D]= ndgrid(Ky_cut,W_cut);
    % K_cut_2D = sqrt(Kx_cut_2D.^2 + Ky_cut_2D.^2);
    K_cut_2D = repmat(K_cut,[size(W_cut) 1])';

    % get power spectra in cut KW
    power_cut_dir = F(Kx_cut_2D,Ky_cut_2D,W_cut_2D);
    power_cut_dir_log = 10*log10(power_cut_dir);  

    % get filtered power spectra in cut KW
    power_cut_dir_masked = F_mask(Kx_cut_2D,Ky_cut_2D,W_cut_2D);
    power_cut_dir_masked_log = 10*log10(power_cut_dir_masked); 

    % get dispersion relation for K_cut
    w_DS_pos = retrieve_w_from_dispersion_relation(Kx_cut,Ky_cut,Ux,Uy,water_depth,1);
    w_DS_neg = retrieve_w_from_dispersion_relation(Kx_cut,Ky_cut,Ux,Uy,water_depth,-1);
    w_DS_merged = [w_DS_pos w_DS_neg];
    K_cut_merged = [K_cut K_cut];

    %%%%%%%%
    % plot
    %%%%%%%%

    % plot W average spectrum
    axes(ax(1));
    hold off;
    pcolor(Spectrum.Kx_3D(:,:,1),Spectrum.Ky_3D(:,:,1),power_Spectrum_2D_DR_log);
    shading flat;
    axis xy equal tight;
    hold on;
    plot(Kx_cut,Ky_cut,'r','Linewidth',2);
    xlabel('Wave number [rad/m]');
    ylabel('Wave number [rad/m]');
    legend('power AV',[num2str(w_cut_dir) ' [deg]' ]);
    title('W average');
    text(Kx_cut(1),Ky_cut(1),'   (-)','Color',[1 0 0]);
    text(Kx_cut(end),Ky_cut(end),'   (+)','Color',[1 0 0]);


    % plot K W cut
    axes(ax(2));
    hold off;
    % offset pcolor
    ox = (K_cut(2) - K_cut(1)) / 2;
    oy = (W_cut(2) - W_cut(1)) / 2;

    pcolor(K_cut_2D - ox ,W_cut_2D - oy ,power_cut_dir_log);
    shading flat;
    axis xy tight;
    xlabel('Wave number [rad/m]');
    ylabel('Wave frequency [rad/sec]');
    hold on;
    plot(K_cut_merged,w_DS_merged,'.k');
    % plot(K_cut,w_DS_min,'r');
    % plot(K_cut,w_DS_max,'g');
    title({['Dispersion Relation cut ' num2str(w_cut_dir) ' [deg]' ],[' Ux: ' num2str(Ux) ' [ m/s]'],...
          [' Uy: ' num2str(Uy) ' [ m/s]']});

    xlim([min(K_cut_2D(:)) max(K_cut_2D(:))]-ox);
    ylim([min(W_cut_2D(:)) max(W_cut_2D(:))]-oy); 



    % plot K W cut filtered spectra
    axes(ax(4));
    hold off;

    pcolor(K_cut_2D - ox ,W_cut_2D - oy ,power_cut_dir_masked_log);
    shading flat;
    axis xy tight;
    xlabel('Wave number [rad/m]');
    ylabel('Wave frequency [rad/sec]');
    hold on;
    % plot(K_cut,w_DS,'k');
    plot(K_cut_merged,w_DS_merged,'.k');
    title({['Filter width: ' num2str(w_width) ' [rad/sec]'],[' Ux: ' num2str(Ux) ' [ m/s]'],...
          [' Uy: ' num2str(Uy) ' [ m/s]']});

    xlim([min(K_cut_2D(:)) max(K_cut_2D(:))]-ox);
    ylim([min(W_cut_2D(:)) max(W_cut_2D(:))]-oy);
    
    
    % refresh annotation
    h.UserData.h_ann_1.String =  ['filter width ' num2str(w_width) ' [rad/s]' ];
    h.UserData.h_ann_2.String =  ['W cut direction: ' num2str(w_cut_dir) ' [deg]' ];
    h.UserData.h_ann_3.String =  ['Uy: ' num2str(Uy) ' [m/s]' ];
    h.UserData.h_ann_4.String =  ['Ux: ' num2str(Ux) ' [m/s]' ];
    
    disp('Plot uptdated');
    
end



