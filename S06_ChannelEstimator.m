%---------------------------------------------------------------------------
%  S06_ChannelEstimator
%  Estimate the channel response H at each mid-amble via regularized
%  (Tikhonov) least squares in the frequency domain, then window the CIR to
%  the measured delay-spread support (Wpre/Wpost) to discard noise outside
%  the true support.
%
%  Input : y       - Doppler-compensated complex baseband (from S05_DopplerEstimator)
%          win_smp - per-anchor sample marker (S02_ReferenceBuilder)
%          XPm     - mid-amble reference spectrum (S02_ReferenceBuilder)
%          del     - Tikhonov regularizer (S02_ReferenceBuilder)
%          P       - system parameter struct (Lp, L_block, Wpre, Wpost)
%  Output: Hb      - channel estimate grid at each anchor [n_mid x L_block]
%---------------------------------------------------------------------------
function Hb = S06_ChannelEstimator(y,win_smp,XPm,del,P)

n = numel(win_smp);
Hb = complex(zeros(n,P.L_block));
for m = 1:n
    % Regularized spectral division at anchor m
    seg = y(win_smp(m):win_smp(m)+P.Lp-1);
    H = fft(seg).*conj(XPm(m,:))./(abs(XPm(m,:)).^2+del);
    h = ifft(H);
    % Window the CIR to the measured delay-spread support, then FFT back to L_block
    hp = complex(zeros(1,P.L_block));
    hp(1:P.Wpost) = h(1:P.Wpost);
    hp(P.L_block-P.Wpre+1:P.L_block) = h(P.Lp-P.Wpre+1:P.Lp);
    Hb(m,:) = fft(hp);
end
