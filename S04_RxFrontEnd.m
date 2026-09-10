%---------------------------------------------------------------------------
%  S04_RxFrontEnd
%  Full-band filtering, cut the frame to length starting from the known
%  SYNC onset (frame timing assumed known), split into the f0/f1 branches,
%  and IQ downconvert to baseband (LPF 400 Hz).
%
%  Input : rx  - received signal (after channel + AWGN)
%          F   - frame layout struct (S01_FrameBuilder)
%          P   - system parameter struct (S00_SystemParams)
%  Output: y0,y1 - complex baseband of the f0 and f1 branches
%---------------------------------------------------------------------------
function [y0,y1] = S04_RxFrontEnd(rx,F,P)

% Full-band filter around [f0-df, f1+df]
y = filtfilt(P.b_all,P.a_all,rx);

% Cut the frame to exact length from the known SYNC onset
iy = P.idx_sim + P.N*(numel(P.sync)+P.z+numel(P.preamble)+F.n_tx_bits+P.z) - 1;
if iy > numel(y)
    yx = [y(P.idx_sim:end), zeros(1,iy-numel(y))];
else
    yx = y(P.idx_sim:iy);
end

% Split into the f0/f1 branches (BW=800Hz), then IQ downconvert to baseband
q0 = filtfilt(P.b0f,P.a0f,yx);
q1 = filtfilt(P.b1f,P.a1f,yx);
t = (0:numel(yx)-1)/P.fs;
y0 = filtfilt(P.b_lpf,P.a_lpf,q0.*cos(2*pi*(P.fc-P.df)*t)) + ...
     1j*filtfilt(P.b_lpf,P.a_lpf,q0.*sin(2*pi*(P.fc-P.df)*t));
y1 = filtfilt(P.b_lpf,P.a_lpf,q1.*cos(2*pi*(P.fc+P.df)*t)) + ...
     1j*filtfilt(P.b_lpf,P.a_lpf,q1.*sin(2*pi*(P.fc+P.df)*t));
