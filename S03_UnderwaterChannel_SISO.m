%---------------------------------------------------------------------------
%  S03_UnderwaterChannel_SISO
%  DoVietHa shallow underwater acoustic channel model: sum of reflected
%  paths, each delayed, randomly phase-shifted, and Doppler-shifted by a
%  pure frequency offset (no time-scaling).
%
%  Input : sa   - transmit signal as an analytic signal
%          f_n  - per-path Doppler shift (Hz), already scaled to the target fD
%          th   - per-path random phase (rad), one value per realization
%          Nf   - frame length (samples)
%          P    - system parameter struct (c_norm, k_delays, fs)
%  Output: rx   - received signal before AWGN is added
%---------------------------------------------------------------------------
function rx = S03_UnderwaterChannel_SISO(sa,f_n,th,Nf,P)

n_vec = 0:Nf-1;
rx = zeros(1,Nf);
% Sum over paths: delay + Doppler shift + random phase + normalized amplitude
for p = 1:numel(P.c_norm)
    kp = P.k_delays(p);
    if kp > 0
        sdl = [zeros(1,kp), sa(1:Nf-kp)];
    else
        sdl = sa;
    end
    rx = rx + P.c_norm(p)*real(exp(1j*(2*pi*f_n(p)*n_vec/P.fs+th(p))).*sdl);
end
