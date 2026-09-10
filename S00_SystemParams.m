%---------------------------------------------------------------------------
%  S00_SystemParams
%  All locked design parameters of the system: two-tone BFSK, branch filter
%  bank, mid-amble frame structure, DoVietHa channel, and Monte-Carlo run
%  schedule. Independent of message content -> call once at the top of Main.m.
%
%  Input : MODE  - 'smoke' (quick test) or 'full' (official results)
%  Output: P     - struct holding every system parameter
%---------------------------------------------------------------------------
function P = S00_SystemParams(MODE)

% BFSK signal parameters (paper HUST + branch filter per the original fsk_rxx.m)
P.fs = 192000; P.fc = 12400; P.df = 400; P.br = 500;
P.N  = round(P.fs/P.br);      % samples per bit
P.z  = 50;                    % guard interval (bits)
P.BW = 800;                   % branch filter bandwidth (Hz)

% Mid-amble structure: 7-bit m-sequence + cyclic prefix/suffix
P.pn7 = [1 1 1 0 0 1 0];
P.L_cp = 2; P.L_pn = 7; P.L_cs = 6;
P.mid_bits = [P.pn7(end-P.L_cp+1:end), P.pn7, P.pn7(1:P.L_cs)];
P.L_mid = numel(P.mid_bits);
P.K_data = 22;                % data bits between two mid-ambles
P.K_sub  = 20;                % FDE equalization block size
P.OLS_MARGIN = 1.25;          % margin factor for Overlap-Save overlap
P.ETA_K = 0.125;              % weight of the channel-uncertainty term in the MMSE regularizer

% Preamble (anti-shock) and SYNC (frame acquisition + Doppler estimation)
P.preamble = repmat([1 0],1,32);
P.sync = [1 1 1 1 1 1 1 1 0 1 1 1 0 0 0 1 0 0 1 1 1 0 1 0 0 1 0 1 1 0 0 1 ...
    1 1 1 1 0 1 0 0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 1 1 1 0 1 1 1 1 0 0 ...
    1 0 0 0 0 1 0 1 1 1 0 1 1 0 1 0 1 0 0 0 1 1 1 0 0 0 0 1 1 1 0 1 ...
    1 1 0 0 0 1 0 1 0 0 1 1 1 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 1 0 1 1 ...
    1 0 1 0 1 1 0 1 0 0 0 1 1 1 0 0 0 0 1 1 0 1 1 0 1 1 0 0 0 1 0 1 ...
    0 0 1 1 1 1 1 1 0 1 1 0 0 1 0 0 0 0 1 1 1 0 1 1 0 1 0 0 0 1 1 0 ...
    1 0 1 0 1 1 0 0 1 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 1 1 0 1 0 1 0 1 ...
    0 1 0 0 1 1 0 1 1 1 0 0 0 1 1 0 1 1 0 0 1 1 0];
P.n_sync_smp = numel(P.sync)*P.N;
P.n_skip_smp = P.N*(numel(P.sync)+P.z+numel(P.preamble));   % samples skipped before the payload
P.idx_sim    = P.N*numel(P.preamble)+1;                     % sample where SYNC starts (frame timing assumed known)

% SYNC sub-blocks used for Doppler estimation (unambiguous range +-8.33 Hz)
P.L_sb = 30;
P.n_sb = floor(numel(P.sync)/P.L_sb);
P.dt_sb = P.L_sb*P.N/P.fs;

% Filters: full band, the two branches f0/f1, and the baseband LPF
[P.b_lpf,P.a_lpf] = butter(4,P.df/(P.fs/2),'low');
[P.b_all,P.a_all] = butter(4,[(P.fc-P.df-P.br)/(P.fs/2),(P.fc+P.df+P.br)/(P.fs/2)],'bandpass');
[P.b0f,P.a0f] = butter(4,[(P.fc-P.df-P.BW/2)/(P.fs/2),(P.fc-P.df+P.BW/2)/(P.fs/2)],'bandpass');
[P.b1f,P.a1f] = butter(4,[(P.fc+P.df-P.BW/2)/(P.fs/2),(P.fc+P.df+P.BW/2)/(P.fs/2)],'bandpass');

% DoVietHa channel: 30 reflected paths, delay + native Doppler + amplitude
PARAMS_FILE = 'c:\Users\Acer\WICOM\SISO\UW-Siso-DoVietHa2017\24_8_26\channel_params_output.mat';
cp = load(PARAMS_FILE);
P.c_n = cp.c_n_proposed; P.tau_n = cp.tau_n_proposed;
P.L_ch = cp.L; P.f_n0 = cp.f_n_proposed;
tau_rel = P.tau_n - min(P.tau_n);
P.k_delays = round(tau_rel*P.fs);
P.c_norm = P.c_n/sqrt(sum(P.c_n.^2));
P.N_cir = max(P.k_delays)+1;

% Mid-amble FFT length and the measured -40dB impulse-response support
P.Lp = P.L_pn*P.N;
Ltest = 80*P.N; tt = (0:Ltest-1)/P.fs; b_st = 40*P.N+1; b_en = 41*P.N;
imp = zeros(1,Ltest); imp(b_st:b_en) = 1;
sg = cos(2*pi*(P.fc-P.df)*tt).*imp;
yy = filtfilt(P.b_all,P.a_all,sg); qq = filtfilt(P.b0f,P.a0f,yy);
uu = filtfilt(P.b_lpf,P.a_lpf,qq.*cos(2*pi*(P.fc-P.df)*tt)) + ...
     1j*filtfilt(P.b_lpf,P.a_lpf,qq.*sin(2*pi*(P.fc-P.df)*tt));
ee = abs(uu).^2; ee = ee/max(ee); ixx = find(ee>1e-4);
P.Wpre = max(0,b_st-ixx(1));
P.Wpost = max(0,ixx(end)-b_en)+P.N_cir;
P.Npre = ceil(P.OLS_MARGIN*P.Wpre);
P.Npost = ceil(P.OLS_MARGIN*P.Wpost);
if P.Wpre+P.Wpost > P.Lp
    error('SystemParams:support','Lp < measured support: design infeasible.');
end

% Overlap-Save FFT block size and in-band bin index
P.L_block = 2^nextpow2(P.K_sub*P.N+P.Npre+P.Npost);
kb_b = round(P.df*P.L_block/P.fs);
P.inband = [1:kb_b+1, P.L_block-kb_b+1:P.L_block];

% Anchor rate and the design limit from the rule R_anchor > 3.35*fD
P.R_anchor = P.fs/((P.L_mid+P.K_data)*P.N);
P.fD_design_limit = P.R_anchor/3.35;

% Monte-Carlo schedule: realization budget allocated adaptively per fD
% (low fD needs more realizations since its BER is lower and harder to resolve)
if strcmpi(MODE,'full')
    P.fD_vec = [0 1 2 3 4 6 8 16];
    P.Number_Relz = [140 140 120 80 70 30 30 30];
    P.snr_vec = 0:5:30;
else
    P.fD_vec = [0 4 16];
    P.Number_Relz = [5 5 5];
    P.snr_vec = [0 10 20 30];
end
