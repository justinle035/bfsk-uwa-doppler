%---------------------------------------------------------------------------
%  S02_ReferenceBuilder
%  Build the known reference signals the receiver needs to correlate
%  against: the SYNC reference IQ (used for Doppler estimation) and the
%  per-anchor mid-amble reference spectrum (used for channel estimation).
%
%  Input : P  - system parameter struct (S00_SystemParams)
%          F  - frame layout struct (S01_FrameBuilder)
%  Output: Rf - reference struct: xp0s/xp1s, win_smp/win_ctr, XP0m/XP1m, del0/del1
%---------------------------------------------------------------------------
function Rf = S02_ReferenceBuilder(P,F)

% SYNC reference IQ (used as the reference for the Doppler estimator)
t_sr = (0:P.n_sync_smp-1)/P.fs;
xs0 = cos(2*pi*(P.fc-P.df)*t_sr).*linecoding(P.sync,P.N,0);
xs1 = cos(2*pi*(P.fc+P.df)*t_sr).*linecoding(P.sync,P.N,1);
Rf.xp0s = filtfilt(P.b_lpf,P.a_lpf,xs0.*cos(2*pi*(P.fc-P.df)*t_sr)) + ...
          1j*filtfilt(P.b_lpf,P.a_lpf,xs0.*sin(2*pi*(P.fc-P.df)*t_sr));
Rf.xp1s = filtfilt(P.b_lpf,P.a_lpf,xs1.*cos(2*pi*(P.fc+P.df)*t_sr)) + ...
          1j*filtfilt(P.b_lpf,P.a_lpf,xs1.*sin(2*pi*(P.fc+P.df)*t_sr));

% Sample marker of each mid-amble (after the cyclic prefix L_cp) in the received frame
mid_smp = P.n_skip_smp + (F.mid_bit_pos-1)*P.N + 1;
Rf.win_smp = mid_smp + P.L_cp*P.N;
Rf.win_ctr = Rf.win_smp + P.Lp/2;

% Mid-amble reference spectrum at each anchor (m-sequence repeated 3x for stable filtering)
Rf.XP0m = zeros(F.n_mid,P.Lp); Rf.XP1m = zeros(F.n_mid,P.Lp);
b3_0 = linecoding(repmat(P.pn7,1,3),P.N,0);
b3_1 = linecoding(repmat(P.pn7,1,3),P.N,1);
for m = 1:F.n_mid
    t0 = (Rf.win_smp(m)-1-P.Lp)/P.fs;
    t3 = t0+(0:3*P.Lp-1)/P.fs;
    r0 = cos(2*pi*(P.fc-P.df)*t3).*b3_0;
    r1 = cos(2*pi*(P.fc+P.df)*t3).*b3_1;
    i0 = filtfilt(P.b_lpf,P.a_lpf,r0.*cos(2*pi*(P.fc-P.df)*t3)) + ...
         1j*filtfilt(P.b_lpf,P.a_lpf,r0.*sin(2*pi*(P.fc-P.df)*t3));
    i1 = filtfilt(P.b_lpf,P.a_lpf,r1.*cos(2*pi*(P.fc+P.df)*t3)) + ...
         1j*filtfilt(P.b_lpf,P.a_lpf,r1.*sin(2*pi*(P.fc+P.df)*t3));
    Rf.XP0m(m,:) = fft(i0(P.Lp+1:2*P.Lp));
    Rf.XP1m(m,:) = fft(i1(P.Lp+1:2*P.Lp));
end

% Tikhonov regularizer for spectral division at the mid-amble
kb_p = round(P.df*P.Lp/P.fs);
inb_p = [1:kb_p+1, P.Lp-kb_p+1:P.Lp];
Rf.del0 = 1e-3*mean(abs(Rf.XP0m(1,inb_p)).^2);
Rf.del1 = 1e-3*mean(abs(Rf.XP1m(1,inb_p)).^2);
