function [h1,h2] = plot_STCFIT_multi_K(STCFIT_multi_K,window_num_vec2plot, plot_video_flag)
% plot K dependence in fit


if exist('plot_video_flag','var') == 0 || isempty(plot_video_flag)
    plot_video_flag = 0;
end

% get number of window 
window_num_vec = [STCFIT_multi_K(:).out_fit.window_used];

% get index to plot
ind = find(window_num_vec == window_num_vec2plot,1);

if ~isempty(ind)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  plot depth fit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

	h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    
%     ax(1) = subplot(1,3,1);
%     scatter3(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep, ...
%             STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep,...
%             STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV,...
%             50,STCFIT_multi_K.out_fit(1).SG_fit.SNR_density_max_Kstep,'filled');
%     
%     xlabel('Ux [m/s]');  
%     ylabel('Uy [m/s]'); 
%     zlabel('K [wavenumber]'); 
%     % axis xy equal tight;
%     
%     hc = colorbar;    
%     ylabel(hc,'SNR');
%     
%     xlim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:))]);
%     ylim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:))]);
%     
%     ax(2) = subplot(1,3,2);
%     
%     wavelength = 2*pi./STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV;
%     scatter3(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep, ...
%             STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep,...
%             wavelength,...
%             50,STCFIT_multi_K.out_fit(1).SG_fit.SNR_density_max_Kstep,'filled');
%     
%     xlabel('Ux [m/s]');  
%     ylabel('Uy [m/s]'); 
%     zlabel('Wavelength [m]'); 
%     axis xy equal tight;
%     
%     hc = colorbar;    
%     ylabel(hc,'SNR');
%     
%     xlim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:))]);
%     ylim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:))]);
    
    nind = 1;
    ax(nind) = subplot(1,3,nind);
    
    wavelength = 2*pi./STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV;
    % depth = 0.078 * wavelength; (7.8%)
    depth = 0.078 * wavelength;
    scatter3(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep, ...
            STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep,...
            depth,...
            50,STCFIT_multi_K.out_fit(1).SG_fit.SNR_density_max_Kstep,'filled');
    
    xlabel('Ux [m/s]');  
    ylabel('Uy [m/s]'); 
    zlabel('Depth [m]'); 
    % axis xy equal tight;
    
    hc = colorbar;    
    ylabel(hc,'SNR');
    
    xlim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:))]);
    ylim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:))]);
    
    
    nind = 2;
    ax(nind) = subplot(1,3,nind);
    
    wavelength = 2*pi./STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV;
    % depth = 0.078 * wavelength; (7.8%)
    depth = 0.078 * wavelength;
    scatter(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep, ...
            depth,...
            50,STCFIT_multi_K.out_fit(1).SG_fit.SNR_density_max_Kstep,'filled');
    
    xlabel('Ux [m/s]');  
    ylabel('Depth [m]'); 
    % axis xy equal tight;
    
    hc = colorbar;    
    ylabel(hc,'SNR');
    
    xlim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:))]);
    ylim([min(depth(:)) max(depth(:))]);
    
    
    nind = 3;
    ax(nind) = subplot(1,3,nind);
    
    wavelength = 2*pi./STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV;
    % depth = 0.078 * wavelength; (7.8%)
    depth = 0.078 * wavelength;
    scatter(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep, ...
            depth,...
            50,STCFIT_multi_K.out_fit(1).SG_fit.SNR_density_max_Kstep,'filled');
    
    xlabel('Uy [m/s]');  
    ylabel('Depth [m]'); 
    % axis xy equal tight;
    
    hc = colorbar;    
    ylabel(hc,'SNR');
    
    xlim([min(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:))]);
    ylim([min(depth(:)) max(depth(:))]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  plot SNR first Guess
    %%%%%%%%%%%%%%%%%%%%%%%%%% 
    h2 = figure('units','normalized','outerposition',[0 0 1 1]);

    offset_x_FG =  (STCFIT_multi_K.out_fit(ind).FG_fit.Ux_2D(2,1) - STCFIT_multi_K.out_fit(ind).FG_fit.Ux_2D(1,1))/2;
    offset_y_FG =  (STCFIT_multi_K.out_fit(ind).FG_fit.Uy_2D(1,2) - STCFIT_multi_K.out_fit(ind).FG_fit.Uy_2D(1,1))/2;
    
    xlimits_vec = [min(STCFIT_multi_K.out_fit(ind).FG_fit.Ux_2D(:)) max(STCFIT_multi_K.out_fit(ind).FG_fit.Ux_2D(:))] - offset_x_FG;
    ylimits_vec = [min(STCFIT_multi_K.out_fit(ind).FG_fit.Uy_2D(:)) max(STCFIT_multi_K.out_fit(ind).FG_fit.Uy_2D(:))] - offset_y_FG;
    
    pcolor(STCFIT_multi_K.out_fit(ind).FG_fit.Ux_2D - offset_x_FG , ...
           STCFIT_multi_K.out_fit(ind).FG_fit.Uy_2D - offset_y_FG , ...
           STCFIT_multi_K.out_fit(ind).FG_fit.SNR_density_2D);
       
    shading flat;
    axis xy equal tight;
    hold on;
    plot(STCFIT_multi_K.out_fit(ind).FG_fit.Ux_fit, STCFIT_multi_K.out_fit(ind).FG_fit.Uy_fit,'+k','Linewidth',2);
    % plot(Ux, Uy,'+k','Linewidth',2);
    xlabel('Ux [m/s]');
    ylabel('Uy [m/s]');
    cb3 = colorbar;
    ylabel(cb3,'SNR density')
    
    xlim(xlimits_vec);
    ylim(ylimits_vec);
    
    str_U_FG = ['Ux = ' num2str(STCFIT_multi_K.out_fit(ind).FG_fit.Ux_fit) ' [m/s]' ...
             '   Uy = ' num2str(STCFIT_multi_K.out_fit(ind).FG_fit.Uy_fit) ' [m/s]'];
    title({['First guess: ' str_U_FG ],['SNR density = ' num2str(STCFIT_multi_K.out_fit(ind).FG_fit.SNR_density_max)]});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot video
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plot_video_flag == 1
        
        n_k_fits = size(STCFIT_multi_K.out_fit(ind).SG_fit.signal_3D,3);
        
        h3 = figure('units','normalized','outerposition',[0 0 1 1]);
        
        offset_x_SG =  (STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(2,1) - STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(1,1))/2;
        offset_y_SG =  (STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(1,2) - STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(1,1))/2;

        xlimits_vec = [min(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D(:))] - offset_x_SG;
        ylimits_vec = [min(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:)) max(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D(:))] - offset_y_SG;
        
        % axis_color = prctile(10*log10(STCFIT_multi_K.out_fit(ind).SG_fit.SNR_density_3D(:)),[1 99]);
        % caxis(axis_color);
        
        for i1 = 1:n_k_fits

            pcolor(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_2D - offset_x_SG ,...
                   STCFIT_multi_K.out_fit(ind).SG_fit.Uy_2D - offset_y_SG ,...
                   STCFIT_multi_K.out_fit(ind).SG_fit.SNR_density_3D(:,:,i1));
               
            shading flat;
            axis xy equal tight;
            hold on;
            plot(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep(i1), ...
                 STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep(i1),'+r','Linewidth',2);
            xlabel('Ux [m/s]');
            ylabel('Uy [m/s]');
            cb3 = colorbar;
            ylabel(cb3,'SNR density')

            xlim(xlimits_vec);
            ylim(ylimits_vec);

            str_U_SG = ['Ux = ' num2str(STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep(i1)) ' [m/s]' ...
                     '   Uy = ' num2str(STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep(i1)) ' [m/s]'];
                 
            title({['Second guess: ' str_U_SG ],...
                   ['SNR density = ' num2str(STCFIT_multi_K.out_fit(ind).SG_fit.SNR_density_max_Kstep(i1))],...
                   ['K = ' num2str(STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV(i1)) ...
                   '  Wavelength = '  num2str(2*pi/STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV(i1)) ' m']});
            
            tmp_2D =    STCFIT_multi_K.out_fit(ind).SG_fit.SNR_density_3D(:,:,i1);
            % axis_color = prctile(tmp_2D(:),[1 99]);   
            axis_color = [min(tmp_2D(:)) max(tmp_2D(:))]; 
            
            caxis(axis_color);   
            pause(0.5);   
            % pause;
        end
        
        
    end
    
    
else
   disp('window_num_vec2plot not found in STCFIT_multi_K structure.'); 
   disp(['fitted windows: ' num2str(window_num_vec)]); 
   h = [];
end




end

