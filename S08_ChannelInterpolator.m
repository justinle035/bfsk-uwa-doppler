%---------------------------------------------------------------------------
%  S08_ChannelInterpolator
%  Linearly interpolate the estimated channel between two consecutive
%  anchors to an arbitrary time tc (the center of the data block being
%  equalized).
%
%  Input : Hb  - channel estimate grid at the anchors (S06_ChannelEstimator)
%          ctr - per-anchor time marker in samples (S02_ReferenceBuilder)
%          tc  - target time to interpolate to
%  Output: Hi  - interpolated channel at time tc
%---------------------------------------------------------------------------
function Hi = S08_ChannelInterpolator(Hb,ctr,tc)

n = size(Hb,1);
% Out of range: hold the nearest anchor value
if tc <= ctr(1), Hi = Hb(1,:); return; end
if tc >= ctr(n), Hi = Hb(n,:); return; end
% Linear interpolation between the two anchors surrounding tc
m = find(ctr<=tc,1,'last');
lam = (tc-ctr(m))/(ctr(m+1)-ctr(m));
Hi = (1-lam)*Hb(m,:) + lam*Hb(m+1,:);
