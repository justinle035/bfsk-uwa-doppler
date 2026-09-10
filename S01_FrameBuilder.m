%---------------------------------------------------------------------------
%  S01_FrameBuilder
%  Interleave mid-ambles into the data, assemble preamble+SYNC+guard+payload
%  into the two-tone BFSK transmit signal, and compute the sample markers
%  the receiver will need.
%
%  Input : payload - message bit vector to send
%          P        - system parameter struct (S00_SystemParams)
%  Output: sa       - transmit signal as an analytic (Hilbert) signal
%          Nf       - frame length (samples)
%          F        - frame layout struct: mid-amble/data positions, FDE blocks
%---------------------------------------------------------------------------
function [sa,Nf,F] = S01_FrameBuilder(payload,P)

F.n_bits = numel(payload);

% Interleave mid-ambles (anchors) and data: K_data bits between two anchors
M_seg = ceil(F.n_bits/P.K_data);
payload_tx = []; mid_bit_pos = zeros(1,M_seg+1);
dat_bit_pos = zeros(1,M_seg); dat_len = zeros(1,M_seg);
ptr = 1; dptr = 1;
for m = 1:M_seg
    mid_bit_pos(m) = ptr;
    payload_tx = [payload_tx, P.mid_bits]; %#ok<AGROW>
    ptr = ptr+P.L_mid;
    nk = min(P.K_data,F.n_bits-dptr+1);
    dat_bit_pos(m) = ptr; dat_len(m) = nk;
    payload_tx = [payload_tx, payload(dptr:dptr+nk-1)]; %#ok<AGROW>
    ptr = ptr+nk; dptr = dptr+nk;
end
mid_bit_pos(M_seg+1) = ptr;
payload_tx = [payload_tx, P.mid_bits];     % final anchor, closes the frame

F.payload_tx = payload_tx;
F.mid_bit_pos = mid_bit_pos;
F.dat_bit_pos = dat_bit_pos;
F.dat_len = dat_len;
F.M_seg = M_seg;
F.n_mid = M_seg+1;
F.n_tx_bits = numel(payload_tx);

% Assemble the bits into the two-tone BFSK signal (preamble+SYNC+guard+payload_tx+guard)
ze = zeros(1,P.z*P.N);
b0 = [linecoding(P.preamble,P.N,0) linecoding(P.sync,P.N,0) ze ...
      linecoding(P.preamble,P.N,0) linecoding(payload_tx,P.N,0) ze];
b1 = [linecoding(P.preamble,P.N,1) linecoding(P.sync,P.N,1) ze ...
      linecoding(P.preamble,P.N,1) linecoding(payload_tx,P.N,1) ze];
t = (0:numel(b0)-1)/P.fs;
s = cos(2*pi*(P.fc-P.df)*t).*b0 + cos(2*pi*(P.fc+P.df)*t).*b1;
s = s/max(abs(s));
Nf = numel(s);
sa = hilbert(s);

% Split each data segment into K_sub-bit FDE equalization blocks
sub = struct('start',{},'nb',{},'first_bit',{});
gb = 1;
for m = 1:M_seg
    nrem = dat_len(m); off = 0;
    while nrem > 0
        nb = min(P.K_sub,nrem);
        sub(end+1) = struct('start',P.n_skip_smp+(dat_bit_pos(m)-1+off)*P.N+1, ...
            'nb',nb,'first_bit',gb); %#ok<AGROW>
        gb = gb+nb; off = off+nb; nrem = nrem-nb;
    end
end
F.sub = sub;
