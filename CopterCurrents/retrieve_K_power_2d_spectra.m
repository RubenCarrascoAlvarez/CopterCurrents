function [spectra_2d,Spectrum,power_Spectrum_3D_DR,h] = retrieve_K_power_2d_spectra(IMG_SEQ,STCFIT,window_number,plot_flag)



%%%%%%%%%%%%%%%%%
% get 3D spectrum
%%%%%%%%%%%%%%%%%

% cut window_number
IMG_SEQ_Window = double(IMG_SEQ.IMG(STCFIT.Windows.w_corners_dim1(2,window_number):STCFIT.Windows.w_corners_dim1(3,window_number),...
                                    STCFIT.Windows.w_corners_dim2(2,window_number):STCFIT.Windows.w_corners_dim2(3,window_number),:));                            
% get spectrum
Spectrum = retrieve_power_spectrum(IMG_SEQ_Window,IMG_SEQ.dx, IMG_SEQ.dy, IMG_SEQ.dt,STCFIT.fit_param.K_limits, STCFIT.fit_param.W_limits); 

%%%%%%%%%%%
% get Ux Uy
%%%%%%%%%%%

% check if STCFIT has been already fitted

if isempty(STCFIT.out_fit)
    disp(['fitting widow ' num2str(window_number) ' ...']);
    
    % fit Spectrum
    out_fit_window_number =  fit_Spectrum2dispersionRelation(Spectrum,STCFIT.fit_param,STCFIT.Windows.average_depth_win(window_number));
    
    Ux_fit = out_fit_window_number.SG_fit.Ux_fit;
    Uy_fit = out_fit_window_number.SG_fit.Uy_fit;
    

    disp(['fitting finished !!!']);
else
    Ux_fit = STCFIT.out_fit(window_number).SG_fit.Ux_fit;
    Uy_fit = STCFIT.out_fit(window_number).SG_fit.Uy_fit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve 2D K spectra
%%%%%%%%%%%%%%%%%%%%%%%%%

% get dispersion realtion
[DR_3D_mask] = get_Dispersion_Relation_3D_mask(Spectrum, ....
                STCFIT.Windows.average_depth_win(window_number), ...
                Ux_fit, Uy_fit, STCFIT.fit_param.w_width_SG);

% get power average in K
power_Spectrum_3D_DR = Spectrum.power_Spectrum;
power_Spectrum_3D_DR(DR_3D_mask~=1) = NaN;
power_Spectrum_2D_DR = nansum(power_Spectrum_3D_DR,3);

Kx_2D= Spectrum.Kx_3D(:,:,1);
Ky_2D= Spectrum.Ky_3D(:,:,1);

% normalize spectrum
power_Spectrum_K_2d = power_Spectrum_2D_DR/sum(power_Spectrum_2D_DR(:));

spectra_2d = struct('power_Spectrum_K_2d',power_Spectrum_K_2d,'Kx_2D',Kx_2D,'Ky_2D',Ky_2D);

if plot_flag == 1
    
    wavelength_x = 2*pi./Kx_2D;
    wavelength_y = 2*pi./Ky_2D;
    log_power_Spectrum_K_2d = 10*log10(spectra_2d.power_Spectrum_K_2d);
    log_power_Spectrum_K_2d(~isfinite(log_power_Spectrum_K_2d)) = NaN;
    
    color_axis = prctile(log_power_Spectrum_K_2d(~isnan(log_power_Spectrum_K_2d)),[0.1 99.9]);

    h = figure('units','normalized','outerposition',[0 0 1 1]);
    colormap('jet');
    
%     ax(1) = subplot(1,2,1);
    pcolor(Kx_2D,Ky_2D,log_power_Spectrum_K_2d);
    shading flat;
    caxis(color_axis);
    axis xy equal tight;
    xlabel('Wavenumber [rad/m]')
    ylabel('Wavenumber [rad/m]')
    grid on
    hc1 = colorbar;
    ylabel(hc1,'power in logaritmic scale [uncalibrated]');
    
    
    % plot axis
    min_kx = min(Kx_2D(:));
    min_ky = min(Ky_2D(:));
    max_kx = max(Kx_2D(:));
    max_ky = max(Ky_2D(:));
    
    hold on;
    plot([0 0] ,[min_ky max_ky],':k','LineWidth',2.5);
    plot([min_kx max_kx], [0 0] ,':k','LineWidth',2.5);
    
    
%     ax(2) = subplot(1,2,2);
%     
%     pcolor(wavelength_x,wavelength_y,log_power_Spectrum_K_2d);
%     shading flat;
%     axis xy ;
%     xlabel('Wavelength [m]')
%     ylabel('Wavelength [m]')
%     grid on
%     hc2 = colorbar;
%     ylabel(hc2,'power in logaritmic scale [uncalibrated]');
    
    FontSize = 18;
    set(findall(h,'-property','FontSize'),'FontSize',FontSize);
    set(findall(h,'-property','FontWeight'),'FontWeight','bold');
    % set(findall(sch,'-property','MarkerEdgeColor'),'MarkerEdgeColor',[0.1 .1 .1]);
    % set(findall(sch,'-property','MarkerEdgeColor'),'LineWidth',2.5);

else
    h = [];
end

            
            
end

