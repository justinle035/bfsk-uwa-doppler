%---------------------------------------------------------------------------
%  S05_DopplerEstimator
%  Estimate the common Doppler shift from SYNC using a Kay-type lag-1
%  phase-increment estimator over 30-bit sub-blocks. Unambiguous range:
%  +-1/(2*dt_sb) = +-8.33 Hz.
%
%  Input : y0,y1 - received complex baseband (from S04_RxFrontEnd)
%          Rf    - reference struct (S02_ReferenceBuilder), uses xp0s/xp1s
%          P     - system parameter struct (L_sb, N, n_sb, dt_sb)
%  Output: f_hat - estimated Doppler frequency (Hz)
%---------------------------------------------------------------------------
function f_hat = S05_DopplerEstimator(y0,y1,Rf,P)

Lsb = P.L_sb*P.N;
s0 = zeros(1,P.n_sb); s1 = zeros(1,P.n_sb);
% Correlation integral over each 30-bit sub-block
for j = 1:P.n_sb
    a = (j-1)*Lsb+1; b = j*Lsb;
    s0(j) = sum(y0(a:b).*conj(Rf.xp0s(a:b)));
    s1(j) = sum(y1(a:b).*conj(Rf.xp1s(a:b)));
end
% Lag-1 phase increment: no unwrap needed, so it never slips a cycle
acc = sum(s0(2:end).*conj(s0(1:end-1))) + sum(s1(2:end).*conj(s1(1:end-1)));
f_hat = angle(acc)/(2*pi*P.dt_sb);
