%---------------------------------------------------------------------------
%  Main
%  1-CLICK RUN: build the transmit frame, run the Monte-Carlo BER sweep over
%  SNR and Doppler fD, save the results for Plot.m to draw. No other file
%  needs to be run first -- open this file, press Run (F5), and it runs the
%  whole S00..S10 pipeline.
%
%  To change the simulation time: edit MODE ('smoke'|'full') below, or edit
%  fD_vec/snr_vec/Number_Relz directly in S00_SystemParams.m.
%---------------------------------------------------------------------------
clear; close all; clc;

MODE = 'smoke';   % 'smoke' = quick test, a few minutes | 'full' = official results, ~30-50 min

fprintf('=========================================================\n');
fprintf(' BFSK/FDE UWA + Doppler tracking -- BER sweep [%s]\n',upper(MODE));
fprintf('=========================================================\n\n');

SCRIPT_DIR = fileparts(mfilename('fullpath'));
out_dir = fullfile(SCRIPT_DIR,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% Set up the system parameters and the message content to send
P = S00_SystemParams(MODE);
payload = text2bits(['trong moi truong nghien cuu khoa hoc hien dai viec ung dung ' ...
    'cac thuat toan vao xu ly tin hieu so la dieu vo cung can thiet ' ...
    'de toi uu hoa hieu suat truyen thong tin tren cac kenh truyen ' ...
    'nhieu nhu am thanh duoi nuoc.']);

% Build the transmit frame and receiver references (once, shared by the whole sweep)
[sa_tx,Nf_tx,F] = S01_FrameBuilder(payload,P);
Rf = S02_ReferenceBuilder(P,F);

fprintf('[CFG] overhead=%.1f%% | anchor rate=%.2f Hz -> design limit fD~%.2f Hz\n',...
    100*(F.n_tx_bits-F.n_bits)/F.n_tx_bits,P.R_anchor,P.fD_design_limit);
fprintf('[CFG] %d anchors, %d FDE blocks, frame length %.2f s\n\n',F.n_mid,numel(F.sub),Nf_tx/P.fs);

% Monte-Carlo sweep: fD x SNR x realization, accumulate raw bit-error counts
nF = numel(P.fD_vec); nS = numel(P.snr_vec);
cnt = zeros(nF,nS); Ntot = zeros(nF,nS);

t_all = tic;
for fi = 1:nF
    fD = P.fD_vec(fi);
    % Scale the native Doppler spectrum to the target fD under test
    if fD > 0
        f_n = P.f_n0*(fD/max(abs(P.f_n0)));
    else
        f_n = zeros(size(P.f_n0));
    end
    tf = tic;
    rng(42);
    for r = 1:P.Number_Relz(fi)
        % Draw one channel realization (delay + Doppler + random phase per path)
        th = 2*pi*rand(1,P.L_ch);
        rx = S03_UnderwaterChannel_SISO(sa_tx,f_n,th,Nf_tx,P);
        for si = 1:nS
            snr = P.snr_vec(si);
            % Add AWGN and run the front end (filtering + IQ downconversion)
            [y0,y1] = S04_RxFrontEnd(awgn(rx,snr,'measured','dB'),F,P);

            % Estimate Doppler from SYNC and compensate the whole received signal
            f_hat = S05_DopplerEstimator(y0,y1,Rf,P);
            t_de = (0:numel(y0)-1)/P.fs;
            dd = exp(-1j*2*pi*f_hat*t_de);
            y0c = y0.*dd; y1c = y1.*dd;

            % Estimate the channel at each anchor and its blind uncertainty
            H0 = S06_ChannelEstimator(y0c,Rf.win_smp,Rf.XP0m,Rf.del0,P);
            H1 = S06_ChannelEstimator(y1c,Rf.win_smp,Rf.XP1m,Rf.del1,P);
            av = 0.5*(S07_ChannelUncertainty(H0,P)+S07_ChannelUncertainty(H1,P));
            reg = 10^(-snr/10)+P.ETA_K*av;   % MMSE regularizer: noise + channel uncertainty

            % Equalize + decide, accumulate the raw bit-error count
            [~,nb] = S10_PayloadDetector(y0c,y1c,F.sub,payload,F.n_bits,H0,H1,Rf.win_ctr,P,reg);
            cnt(fi,si) = cnt(fi,si)+nb;
            Ntot(fi,si) = Ntot(fi,si)+F.n_bits;
        end
    end
    fprintf('--- fD=%4.1f Hz : ',fD); fprintf('%9.2e',cnt(fi,:)./Ntot(fi,:));
    fprintf('  (%.1f s)\n',toc(tf));
end
fprintf('\n[TIME] total %.1f min\n',toc(t_all)/60);

% Save the results for Plot.m to read back and draw (no re-simulation there)
BER = cnt./Ntot;
fD_vec = P.fD_vec; snr_vec = P.snr_vec; R_anchor = P.R_anchor;
save(fullfile(out_dir,'BER_results.mat'),'cnt','Ntot','BER','fD_vec','snr_vec','R_anchor');
fprintf('\n[SAVED] results/BER_results.mat\n[DONE] Main.\n');
