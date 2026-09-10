%---------------------------------------------------------------------------
%  S10_PayloadDetector
%  Process the whole payload: for each FDE block, interpolate the channel to
%  the block center, equalize with Overlap-Save, differentially combine the
%  two branches, make hard decisions per bit, and compute the BER against
%  the original message.
%
%  Input : y0,y1   - Doppler-compensated complex baseband
%          sub     - list of FDE blocks (S01_FrameBuilder)
%          payload - original message (for BER reference)
%          n_bits  - number of message bits
%          H0,H1   - channel estimate grid at the anchors (S06_ChannelEstimator)
%          ctr     - anchor time markers (S02_ReferenceBuilder)
%          P       - system parameter struct
%          reg     - MMSE regularizer (noise + channel uncertainty)
%  Output: ber,nb  - bit error rate and raw error count (for pooling statistics)
%---------------------------------------------------------------------------
function [ber,nb] = S10_PayloadDetector(y0,y1,sub,payload,n_bits,H0,H1,ctr,P,reg)

ms = round(P.N*0.25); me = round(P.N*0.75);   % sampling window at bit center (avoids edge ISI)
bits = zeros(1,n_bits);
for s = 1:numel(sub)
    st = sub(s).start; nbk = sub(s).nb; bl = nbk*P.N; tc = st+bl/2;

    % Interpolate the channel to the block center, use it for the MMSE weights
    A0 = S08_ChannelInterpolator(H0,ctr,tc); A1 = S08_ChannelInterpolator(H1,ctr,tc);
    W0 = conj(A0)./(abs(A0).^2+reg); W1 = conj(A1)./(abs(A1).^2+reg);

    % Overlap-Save equalization, then differentially combine the two branches
    [e0,e1] = S09_Equalizer(y0,y1,st,bl,W0,W1,P);
    d = real(e1)-real(e0);

    % Hard decision per bit in the block
    for b = 0:nbk-1
        k = sub(s).first_bit+b; if k > n_bits, break; end
        bits(k) = mean(d(b*P.N+ms:b*P.N+me)) >= 0;
    end
end
[nb,ber] = biterr(payload(1:n_bits),bits);
