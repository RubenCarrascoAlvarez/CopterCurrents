function [Ux_fit,Uy_fit,depth,Ux,Uy,SNR_density,px,py] = fit_current_depth_profile(STCFIT_multi_K,window_num2use,SNR_density_thr,wavelenth_interval,polyfit_order)
% fit current profile from STCFIT_multi_K
%
% 
% STCFIT_multi_K
% window_num2use: widow fit to use
% SNR_density_thr: SNR threshold for valid fit [2.5];
% wavelenth_interval: wavelenght interval to use in the fit, default [0 inf]
% polifit_order: order of the polifit. [1];


if exist('SNR_density_thr','var') == 0 || isempty(SNR_density_thr)
    SNR_density_thr = 2.5;
end

if exist('wavelenth_interval','var') == 0 || isempty(wavelenth_interval)
    wavelenth_interval = [0 inf];
end

if exist('polyfit_order','var') == 0 || isempty(polyfit_order)
    polyfit_order = 1;
end

% get number of window 
window_num_vec = [STCFIT_multi_K(:).out_fit.window_used];

% get index to use
ind = find(window_num_vec == window_num2use,1);

% get UX,UY in function of K
Ux = STCFIT_multi_K.out_fit(ind).SG_fit.Ux_fit_Kstep;
Uy = STCFIT_multi_K.out_fit(ind).SG_fit.Uy_fit_Kstep;
SNR_density = STCFIT_multi_K.out_fit(ind).SG_fit.SNR_density_max_Kstep;
K = STCFIT_multi_K.out_fit(ind).SG_fit.K_steps_AV;

% apply SNR_density_thr & wavelenth_interval
K_limits = [2*pi/wavelenth_interval(2) 2*pi/wavelenth_interval(1)];

ind2rem  = SNR_density<SNR_density_thr | K< K_limits(1) | K> K_limits(2);

Ux(ind2rem) = [];
Uy(ind2rem) = [];
SNR_density(ind2rem) = [];
K(ind2rem) = [];

wavelength = 2*pi./K;
% depth = 0.078 * wavelength; (7.8%)
depth = 0.078 * wavelength;


% fit in 3D space Ux,Uy,K
% Ux = a1 + b1*depth + c1*depth^2 
% Uy = a2 + b2*depth + c2*depth^2 

% weighted_fit_flag = 0;
weighted_fit_flag = 1;
n = polyfit_order;

fit_flag_depth = 1;
% fit_flag_depth = 0;

% first guess with polinomial fit

if fit_flag_depth == 1

    if weighted_fit_flag == 0
        px = polyfit(depth,Ux,n);
        py = polyfit(depth,Uy,n);
    else
        addpath('/home/carrasco/Matlab/code/polyfitweighted')
        px = polyfitweighted(depth,Ux,n,SNR_density);
        py = polyfitweighted(depth,Uy,n,SNR_density);
    end

    Ux_fit = polyval(px,depth);
    Uy_fit = polyval(py,depth);
    
else
    if weighted_fit_flag == 0
        px = polyfit(K,Ux,n);
        py = polyfit(K,Uy,n);
    else
        addpath('/home/carrasco/Matlab/code/polyfitweighted')
        px = polyfitweighted(K,Ux,n,SNR_density);
        py = polyfitweighted(K,Uy,n,SNR_density);
    end

    Ux_fit = polyval(px,K);
    Uy_fit = polyval(py,K);    
end

% second guess with nlinfit
beta0 = [px(end:-1:1) py(end:-1:1)];
X = double(cat(2,Ux',Uy',depth'));
y = zeros(size(X,1),1);

opts = statset('nlinfit');
opts.MaxIter = 20000;
opts.UseParallel = true;
opts.Display = 'iter'; % 'final'
opts.FunValCheck = 'on';
opts.TolX = 1.0000e-16;

[beta1,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(X,y,@Error_depth_fit_model_order_n,beta0,opts,'ErrorModel','constant');

% compute model
[Ux_fit,Uy_fit] = compute_beta_coef_depth_fit(beta1,depth);

point_size = 55;

h = figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(2,2,1);
scatter3(Ux,Uy,-depth,point_size,SNR_density,'filled');
hold on;
plot3(Ux_fit,Uy_fit,-depth,'--k');
xlabel('Ux [m/s]');  
ylabel('Uy [m/s]'); 
zlabel('Depth [m]'); 
hc1 = colorbar;
ylabel(hc1,'Signal to Noise Ratio');

ax(2) = subplot(2,2,2);
scatter(sqrt(Ux.^2 + Uy.^2),-depth,point_size,SNR_density,'filled');
hold on;
plot(sqrt(Ux_fit.^2 + Uy_fit.^2),-depth,'--k');
xlabel('Uabs [m/s]');  
ylabel('Depth [m]'); 
hc2 = colorbar;
ylabel(hc2,'Signal to Noise Ratio');

ax(3) = subplot(2,2,3);
scatter(Ux,-depth,point_size,SNR_density,'filled');
hold on;
plot(Ux_fit,-depth,'--r');
xlabel('Ux [m/s]');  
ylabel('Depth [m]');
hc3 = colorbar;
ylabel(hc3,'Signal to Noise Ratio');

ax(4) = subplot(2,2,4);
scatter(Uy,-depth,point_size,SNR_density,'filled');
hold on;
plot(Uy_fit,-depth,'--b');
xlabel('Uy [m/s]');  
ylabel('Depth [m]'); 
hc4 = colorbar;
ylabel(hc4,'Signal to Noise Ratio');

colormap('jet');

for i1 = 1:length(Ux_fit)
    disp([ num2str(i1) ' depth: ' num2str(depth(i1)) ' Ux_fit: ' num2str(Ux_fit(i1)) ' Uy_fit: ' num2str(Uy_fit(i1)) ' Uy_fit: ' num2str(sqrt(Ux_fit(i1).^2 + Uy_fit(i1).^2))]);
end

FontSize = 16;
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize);
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold');
set(findall(gcf,'-property','MarkerEdgeColor'),'MarkerEdgeColor',[0.1 .1 .1]);
set(findall(gcf,'-property','MarkerEdgeColor'),'LineWidth',2.5);

l1 = legend(ax(1),'Individual fit','trend','Location','Best');
l2 = legend(ax(2),'Individual fit','trend','Location','Best');
l3 = legend(ax(3),'Individual fit','trend','Location','Best');
l4 = legend(ax(4),'Individual fit','trend','Location','Best');

l1.FontWeight = 'normal'; l1.FontSize = FontSize-3;
l2.FontWeight = 'normal'; l2.FontSize = FontSize-3;
l3.FontWeight = 'normal'; l3.FontSize = FontSize-3;
l4.FontWeight = 'normal'; l4.FontSize = FontSize-3;

end

