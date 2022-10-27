function [f1,f2,p1,p2] = plot_pedm_results(out,figureNumber,U1,U2)

%This function plots the outputs from the 'find_current_depth_profile'
%function for the PEDM.

%INPUTS:
%out - structure output from 'find_current_depth_profile' function
%windowNumber - which spatial window to plot
%figureNumber - which MATLAB figure to plot the results
%U1 (optional): function handle specifying  the true depth profile as a
%function of depth (negative z-values) in the along-flow direction

%U2 (optional): function handle specifying  the true depth profile as a
%function of depth (negative z-values) in the cross-flow direction


zvect = linspace(out.EDM.lin.Z_eff(1),0,400);
figure(figureNumber);f1 = subplot(1,2,1);p1 = plot(...
    out.global.U1,out.EDM.lin.Z_eff,'x',...
    out.PEDM.U1_fun(zvect),zvect,'-k');
hold on;
[U_vert,Z_vert] = calculate_pedm_current_bounds(out.PEDM.U1_fun,out.PEDM.eps1,out.PEDM.verbose1,zvect,0.1);
fill(U_vert,Z_vert,[0,0,0],'EdgeColor','none','FaceAlpha',0.2);
[U_vert,Z_vert] = calculate_pedm_current_bounds(out.PEDM.U1_fun,out.PEDM.eps1,out.PEDM.verbose1,zvect,0.05);
fill(U_vert,Z_vert,[0,0,0],'EdgeColor','none','FaceAlpha',0.4);hold off;    

ylim([min(zvect),0]);
xlabel('[m/s]');ylabel('Depth [m]');
%set(p1(ind1+2),'LineWidth',4.0);
if exist('U1','var')
    hold on;
    plot(U1(zvect),zvect,'--','LineWidth',2.0);
    hold off;
end

% if exist('U1','var')
%     hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','PEDM','True Profile');
%     
% else
%     hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','PEDM');
% end
% set(hl,'Interpreter','latex');
% title('Along-flow');
title('X-direction');



figure(figureNumber);f2 = subplot(1,2,2);p2 = plot(...
    out.global.U2,out.EDM.lin.Z_eff,'x',...
    out.PEDM.U2_fun(zvect),zvect,'-k');
hold on;
[U_vert,Z_vert] = calculate_pedm_current_bounds(out.PEDM.U2_fun,out.PEDM.eps2,out.PEDM.verbose2,zvect,0.1);
fill(U_vert,Z_vert,[0,0,0],'EdgeColor','none','FaceAlpha',0.2);
[U_vert,Z_vert] = calculate_pedm_current_bounds(out.PEDM.U2_fun,out.PEDM.eps2,out.PEDM.verbose2,zvect,0.05);
fill(U_vert,Z_vert,[0,0,0],'EdgeColor','none','FaceAlpha',0.4);hold off;    


if exist('U2','var')
    hold on;
    plot(U2(zvect),zvect,'--','LineWidth',2.0);
    hold off;
end
xlabel('[m/s]');
ylim([min(zvect),0]);

if exist('U2','var')
    hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','PEDM','True Profile');
    
else
    hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','PEDM');
end
set(hl,'Interpreter','latex');
%title('Cross-flow');
title('Y-direction');


end