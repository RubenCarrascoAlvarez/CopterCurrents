function [f1,f2,p1,p2] = plot_depth_profile_results(out,figureNumber,U1,U2)

%This function plots the outputs from the 'find_current_depth_profile'
%function.

%INPUTS:
%out - structure output from 'find_current_depth_profile' function
%figureNumber - which MATLAB figure to plot the results
%U1 (optional): function handle specifying  the true depth profile as a
%function of depth (negative z-values) in the along-flow direction

%U2 (optional): function handle specifying  the true depth profile as a
%function of depth (negative z-values) in the cross-flow direction

zvect = linspace(out.EDM.lin.Z_eff(1),0,400);
[~,ind1] = min([out.PEDM.eps1,out.EDM.log.eps1,out.EDM.lin.eps1]);
[~,ind2] = min([out.PEDM.eps2,out.EDM.log.eps2,out.EDM.lin.eps2]);

figure(figureNumber);f1 = subplot(1,2,1);p1 = plot(...
    out.global.U1,out.EDM.lin.Z_eff,'x',...
    out.global.U1,out.EDM.log.Z_eff,'d',...
    out.PEDM.U1_fun(zvect),zvect,'-k',...
    out.EDM.log.U1_fun(zvect),zvect,'--',...
    out.EDM.lin.U1_fun(zvect),zvect,'-.','LineWidth',2.0);
ylim([min(zvect),0]);
xlabel('[m/s]');ylabel('Depth [m]');
%set(p1(ind1+2),'LineWidth',4.0);
if exist('U1','var')
    hold on;
    plot(U1(zvect),zvect,'LineWidth',4.0);
    hold off;
end


% title(sprintf(...
%     'Along-Flow:\n\\epsilon_{PEDM} = %.4f m/s;\n \\epsilon_{Log} = %.4f m/s, u^* = %.3f m/s;\n \\epsilon_{Lin} = %.4f m/s',...
%     out.PEDM.eps1,out.EDM.log.eps1,out.EDM.log.uStar1,out.EDM.lin.eps1),'horizontalAlignment','left');

title(sprintf(...
    'X-Direction:\n\\epsilon_{PEDM} = %.4f m/s;\n \\epsilon_{Log} = %.4f m/s, u^* = %.3f m/s;\n \\epsilon_{Lin} = %.4f m/s',...
    out.PEDM.eps1,out.EDM.log.eps1,out.EDM.log.uStar1,out.EDM.lin.eps1),'horizontalAlignment','left');

% if exist('U1','var')
%     hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','$\{u_i,Z_{\mathrm{eff,log}}(k_i)\}$',...
%         'PEDM','Log','Linear','True Profile');
%     
% else
%     hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','$\{u_i,Z_{\mathrm{eff,log}}(k_i)\}$',...
%         'PEDM','Log','Linear');
% end
% set(hl,'Interpreter','latex');


figure(figureNumber);f2 = subplot(1,2,2);p2 = plot(...
    out.global.U2,out.EDM.lin.Z_eff,'x',...
    out.global.U2,out.EDM.log.Z_eff,'d',...
    out.PEDM.U2_fun(zvect),zvect,'-k',...
    out.EDM.log.U2_fun(zvect),zvect,'--',...
    out.EDM.lin.U2_fun(zvect),zvect,'-.','LineWidth',2.0);
ylim([min(zvect),0]);
%xlim(sort([0,1.5*median(out.global.U2)]));
xlabel('[m/s]');%ylabel('Depth [m]');

if exist('U2','var')
    hold on;
    plot(U2(zvect),zvect,'LineWidth',4.0);
    hold off;
end


% title(sprintf(...
%     'Cross-Flow:\n\\epsilon_{PEDM} = %.4f m/s;\n \\epsilon_{Log} = %.4f m/s, u^* = %.3f m/s;\n \\epsilon_{Lin} = %.4f m/s',...
%     out.PEDM.eps2,out.EDM.log.eps2,out.EDM.log.uStar2,out.EDM.lin.eps2),'horizontalAlignment','left');
title(sprintf(...
    'Y-direction:\n\\epsilon_{PEDM} = %.4f m/s;\n \\epsilon_{Log} = %.4f m/s, u^* = %.3f m/s;\n \\epsilon_{Lin} = %.4f m/s',...
    out.PEDM.eps2,out.EDM.log.eps2,out.EDM.log.uStar2,out.EDM.lin.eps2),'horizontalAlignment','left');
if exist('U2','var')
    hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','$\{u_i,Z_{\mathrm{eff,log}}(k_i)\}$',...
        'PEDM','Log','Linear','True Profile');
    
else
    hl = legend('$\{u_i,Z_{\mathrm{eff,lin}}(k_i)\}$','$\{u_i,Z_{\mathrm{eff,log}}(k_i)\}$',...
        'PEDM','Log','Linear');
end
set(hl,'Interpreter','latex');
%set(p2(ind2+2),'LineWidth',4.0);








end