function [spectra_1d,Spectrum,power_Spectrum_3D_DR,h] = retrieve_K_power_1d_spectra(IMG_SEQ,STCFIT,window_number,plot_flag)



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

% normalize 2D K spectra
% power_Spectrum_2D_DR_norm = power_Spectrum_2D_DR ./ sum(power_Spectrum_2D_DR(:));
% figure; pcolor(Spectrum.Kx_3D(:,:,1),Spectrum.Ky_3D(:,:,1),log(power_Spectrum_2D_DR)); shading flat;

%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve 1D K spectra
%%%%%%%%%%%%%%%%%%%%%%%%%
K_2D = sqrt(Spectrum.Kx_3D(:,:,1).^2 + Spectrum.Ky_3D(:,:,1).^2);
Kmax = max(K_2D(isfinite(power_Spectrum_2D_DR)));
dK =  sqrt((Spectrum.Kx_3D(2,1,1) - Spectrum.Kx_3D(1,1,1)).^2 + (Spectrum.Ky_3D(1,2,1) - Spectrum.Ky_3D(1,1,1)).^2);
dK_mid = dK/2;
K_1d = 0:dK:Kmax;
power_Spectrum_K_1d = nan(size(K_1d));

for i1 = 1:length(K_1d)
    
    ind2use = K_2D > K_1d(i1)-dK_mid & K_2D < K_1d(i1)+dK_mid;
    if ~isempty(ind2use)
        power_Spectrum_K_1d(i1) = sum(power_Spectrum_2D_DR(ind2use));
    end
end

% figure;plot(K_1d,power_Spectrum_K_1d)

% normalize spectrum
power_Spectrum_K_1d = power_Spectrum_K_1d/sum(power_Spectrum_K_1d(:));

spectra_1d = struct('power_Spectrum_K_1d',power_Spectrum_K_1d,'K_1d',K_1d);

if plot_flag == 1
    
    wavelength = 2*pi./spectra_1d.K_1d;

    h = figure('units','normalized','outerposition',[0 0 1 1]);
    
    ax(1) = subplot(1,2,1);
    bar(spectra_1d.K_1d,spectra_1d.power_Spectrum_K_1d,'b');
    xlabel('Wavenumber [rad/m]')
    ylabel('Power [Uncalibrated]')
    xlim([0 max(spectra_1d.K_1d)/2]);
    grid on
    
    
    ax(2) = subplot(1,2,2);
    
    
    % plot(wave_period,spectra_1d.power_Spectrum_W_1d,'b','LineWidth',2.5);
    % bar(2*pi./spectra_1d.W_1d,spectra_1d.power_Spectrum_W_1d,'b');
    ah = area(wavelength,spectra_1d.power_Spectrum_K_1d,'LineWidth',2.5);
    hold on
    sch = scatter(wavelength,spectra_1d.power_Spectrum_K_1d,40,[1 0 0],'filled');
    xlabel('Wavelength [m]');
    ylabel('Power [Uncalibrated]');
    xlim([0 wavelength(2)/2]);
    grid on
    
    FontSize = 16;
    set(findall(h,'-property','FontSize'),'FontSize',FontSize);
    set(findall(h,'-property','FontWeight'),'FontWeight','bold');
    set(findall(sch,'-property','MarkerEdgeColor'),'MarkerEdgeColor',[0.1 .1 .1]);
    set(findall(sch,'-property','MarkerEdgeColor'),'LineWidth',2.5);
    
    ah.FaceAlpha = 0.3;
    ah.FaceColor = [0 1 0];
    
else
    h = [];
end

            
            
end

