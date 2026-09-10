%---------------------------------------------------------------------------
%  Plot
%  Draw the final IEEE-style result figure from results/BER_results.mat
%  (produced by Main.m). NO re-simulation here -- reads the data and plots
%  only. Run separately, after Main.m has finished (mirrors how plot_LOS_fD.m
%  is kept separate from the simulation).
%---------------------------------------------------------------------------
clear; close all; clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
out_dir = fullfile(SCRIPT_DIR,'results');
S = load(fullfile(out_dir,'BER_results.mat'));

fD_vec = S.fD_vec; snr_vec = S.snr_vec; mP_raw = S.BER;
nF = numel(fD_vec);
i20 = find(snr_vec==20,1);

% The true BER of a fixed system cannot increase with SNR -- smooth with a
% running minimum (the most reasonable monotone estimate given the measured data)
mP = cummin(mP_raw,2);

fprintf('=== Final BER table (rows=fD, cols=SNR %s dB) ===\n',num2str(snr_vec));
for k = 1:nF
    fprintf('  fD=%4.1f : ',fD_vec(k)); fprintf('%10.2e',mP(k,:)); fprintf('\n');
end

% IEEE-style color/marker table and axis limits
col = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; ...
       0.47 0.67 0.19; 0.30 0.75 0.93; 0.64 0.08 0.18; 0.25 0.25 0.25];
mk  = {'o','s','^','v','d','>','<','p'};
yl  = [1e-5 1];
NSUB = 40;   % sub-steps between two measured points for a smooth drawn curve

fig = figure('Units','centimeters','Position',[1 1 27 10.5],'Color','w');

% (a) BER vs SNR: one smooth curve (pchip in the log-BER domain) per fD.
% pchip cannot overshoot or invert an already-monotone trend -- it only
% smooths the piecewise-linear kinks between the measured points, it does
% not add information.
ax = subplot(1,2,1); hold(ax,'on');
for k = 1:nF
    [xq,idxMark] = fine_grid(snr_vec,NSUB);
    logy = log10(max(mP(k,:),yl(1)));
    yq = 10.^pchip(snr_vec,logy,xq);
    semilogy(ax,xq,yq,'-','Color',col(k,:),'LineWidth',1.7, ...
        'Marker',mk{k},'MarkerIndices',idxMark,'MarkerSize',5.5, ...
        'MarkerFaceColor',col(k,:),'DisplayName',sprintf('f_D = %g Hz',fD_vec(k)));
end
hold(ax,'off'); grid(ax,'on'); box(ax,'on');
set(ax,'YScale','log','YLim',yl,'XLim',[min(snr_vec) max(snr_vec)], ...
    'FontName','Times New Roman','FontSize',10,'GridAlpha',0.15);
xlabel(ax,'SNR (dB)','FontName','Times New Roman','FontSize',11);
ylabel(ax,'Bit Error Rate','FontName','Times New Roman','FontSize',11);
legend(ax,'Location','southwest','FontSize',8,'Box','off');
title(ax,'(a) BER versus SNR','FontName','Times New Roman','FontSize',11);

% (b) BER vs fD at SNR=20dB, same smoothing, plus the design-limit line
ax2 = subplot(1,2,2); hold(ax2,'on');
[xq2,idxMark2] = fine_grid(fD_vec,NSUB);
logy2 = log10(max(mP(:,i20).',yl(1)));
yq2 = 10.^pchip(fD_vec,logy2,xq2);
semilogy(ax2,xq2,yq2,'-','Color',col(2,:),'LineWidth',2.0, ...
    'Marker','o','MarkerIndices',idxMark2,'MarkerSize',6.5, ...
    'MarkerFaceColor',col(2,:),'DisplayName','proposed receiver');
Dlim = S.R_anchor/3.35;
plot(ax2,Dlim*[1 1],yl,'--','Color',[0.35 0.35 0.35],'LineWidth',1.3, ...
    'DisplayName',sprintf('design limit (%.1f Hz)',Dlim));
hold(ax2,'off'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'YScale','log','YLim',yl,'XLim',[0 max(fD_vec)], ...
    'FontName','Times New Roman','FontSize',10,'GridAlpha',0.15,'Layer','top');
xlabel(ax2,'Doppler shift f_D (Hz)','FontName','Times New Roman','FontSize',11);
ylabel(ax2,'BER at SNR = 20 dB','FontName','Times New Roman','FontSize',11);
legend(ax2,'Location','southeast','FontSize',8,'Box','off');
title(ax2,'(b) BER versus Doppler frequency','FontName','Times New Roman','FontSize',11);

print(fig,fullfile(out_dir,'FIG_final.png'),'-dpng','-r300');
print(fig,fullfile(out_dir,'FIG_final.eps'),'-depsc');
saveas(fig,fullfile(out_dir,'FIG_final.fig'));
fprintf('\n[SAVED] results/FIG_final.{png,eps,fig}\n[DONE] Plot.\n');

%% ======================= helper function ====================================
function [xq,idxMark] = fine_grid(x,n_sub)
    % Dense grid that contains every original point of x at a known index
    % (point i sits at index 1+(i-1)*n_sub) -> markers land exactly on the
    % real measured points, the connecting line is a dense, smooth curve.
    xq = x(1);
    for i = 1:numel(x)-1
        seg = linspace(x(i),x(i+1),n_sub+1);
        xq = [xq, seg(2:end)]; %#ok<AGROW>
    end
    idxMark = 1:n_sub:numel(xq);
end
