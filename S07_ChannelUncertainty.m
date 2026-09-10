%---------------------------------------------------------------------------
%  S07_ChannelUncertainty
%  Blind channel-uncertainty measure (no true channel needed): how much the
%  channel changes between two consecutive anchors, relative to its own
%  power. This bounds how well linear interpolation can know the channel in
%  between, so it is used as an extra regularizer term in the MMSE equalizer.
%
%  Input : H      - channel estimate grid at the anchors (S06_ChannelEstimator)
%          P      - system parameter struct, uses P.inband
%  Output: d      - channel uncertainty (unitless power ratio)
%---------------------------------------------------------------------------
function d = S07_ChannelUncertainty(H,P)

A = H(:,P.inband);
if size(A,1) < 2
    d = 0; return;
end
% Relative energy difference between consecutive anchors, median for robustness
num = sum(abs(diff(A,1,1)).^2,2);
den = sum(abs(A(1:end-1,:)).^2,2);
d = median(num./max(den,eps));
