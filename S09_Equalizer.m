%---------------------------------------------------------------------------
%  S09_Equalizer
%  Overlap-Save block frequency-domain equalizer core: take extra samples
%  before/after the data block so the FFT equals a true linear convolution,
%  apply the MMSE weights, IFFT, then crop out the equalized data block.
%
%  Input : y0,y1  - Doppler-compensated complex baseband
%          st,bl  - start sample and length (samples) of the data block
%          W0,W1  - MMSE equalizer weights of the two branches (from H + reg)
%          P      - system parameter struct (Npre, Npost, L_block)
%  Output: e0,e1  - equalized signal of the two branches, length bl
%---------------------------------------------------------------------------
function [e0,e1] = S09_Equalizer(y0,y1,st,bl,W0,W1,P)

n0 = numel(y0); en = min(st+bl-1,n0);
b0 = y0(st:en); b1 = y1(st:en); L = numel(b0);

% Take the overlap (before) and lookahead (after) samples around the block
ov = min(P.Npre,st-1);
o0 = y0(st-ov:st-1); o1 = y1(st-ov:st-1);
ps = en+1; la = max(0,min(P.Npost,n0-ps+1));
a0 = y0(ps:ps+la-1); a1 = y1(ps:ps+la-1);
pre = zeros(1,P.Npre-ov); post = zeros(1,P.L_block-P.Npre-L-la);

% FFT -> apply MMSE weights -> IFFT, then crop out the data block
x0 = ifft(fft([pre,o0,b0,a0,post]).*W0);
x1 = ifft(fft([pre,o1,b1,a1,post]).*W1);
v = P.Npre+1;
e0 = zeros(1,bl); e1 = zeros(1,bl);
e0(1:L) = x0(v:v+L-1); e1(1:L) = x1(v:v+L-1);
