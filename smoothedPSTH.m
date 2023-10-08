function [tbins, sm_psth] = smoothedPSTH(spikeTimes, binSize, kernelWindow, alpha,tau)
% This function takes spike time as input and returns time bins and smoothed PSTH 
% [tbins, sm_psth] = smoothedPSTH(spikeTimes, binSize, tWindow, tau, alpha)
% © @Md Rakibul Mowla, Ph.D./University of Iowa
% spikeTimes --> input in seconds 
% Default parameters
if nargin < 2 || isempty(binSize); binSize = 0.05; end % in seconds
if nargin < 3 || isempty(kernelWindow); kernelWindow = 2; end % in seconds 
if nargin < 4 || isempty(tau); tau = 0.1; end 
if nargin < 5 || isempty(alpha); alpha = 20; end 
% spike counts
[Nhist,edges] = histcounts(1000*spikeTimes,'BinWidth',1000*binSize);
sprate = Nhist./binSize;
tbins = edges(1:end-1)./1000;

% causal kernel (alpha function)
kernelT = -kernelWindow/2 : binSize : kernelWindow/2; % kernel length default 2s
kernel = (kernelT.*(alpha^2).*tau).*exp(-alpha.*tau.*kernelT).* (kernelT > 0);
kernel = kernel / sum(kernel); % normalize kernel so its integral is 1

% Convolve firing rates with kernel to obtain smoothed PSTH
sm_psth = conv(sprate, kernel, 'same');

end
