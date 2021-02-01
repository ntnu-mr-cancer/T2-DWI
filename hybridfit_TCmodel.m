function [fitresult, gof] = hybridfit_TCmodel(b, TE, signal, T2_slow_range, T2_fast_range)


% HYBRIDFIT_TCMODEL(B, TE, SIGNAL, T2_SLOW_RANGE, T2_FAST_RANGE)
% Creates a fit of signal fractions for combined T2-DWI,
% using a two-component model with a slow and a fast component.
%
% Model:
% SI/SI_0 = SF_slow exp(-TE/T2_slow) exp(-b*ADC_slow) + SF_fast exp(-TE/T2_fast) exp(-b*ADC_fast)
%         
% and:
% SF_slow + SF_fast = 1 
%
% Data for hybrid fit:
%   b: vector of b-values [b1 b2]
%	TE: vector of TE-values [TE1 TE2]
%	signal: matrix with pixel values [TE1_b1 TE1_b2; TE2_b1 TE2_b2]
%   T2_slow_range: a vector containing the lower and upper limits of
%                   T2_slow [T2_slow_lower T2_slow_upper]
%   T2_fast_range: a vector containing the lower and upper limits of
%                   T2_fast [T2_fast_lower T2_fast_upper]
% Output:
%	fitresult: a fit object representing the fit
%	gof: structure with goodness-of fit info


[bData, teData, signalData] = prepareSurfaceData(b, TE, signal);

% Specifying fixed ADC values for the components (values from literature)
% ADC_slow = 0.3; % µm²/ms
% ADC_fast = 2.6; % µm²/ms

% Setting up fittype and options
ft = fittype('SI_0*((1-SF_fast)*exp(-TE/T2_slow)*exp(-b*0.3)+SF_fast*exp(-TE/T2_fast)*exp(-b*2.6))','independent',{'b','TE'},'dependent','SI','coefficients',{'SF_fast','T2_slow','T2_fast','SI_0'});
opts = fitoptions('Method','NonLinearLeastSquares');
opts.Display = 'Off';
opts.Lower = [0 T2_slow_range(1) T2_fast_range(1) 0]; % SF_fast, T2_slow, T2_fast, SI_0
opts.MaxFunEvals = 10000;
opts.MaxIter = 1000;
opts.StartPoint = [0.5 (T2_slow_range(1)+T2_slow_range(2))/2 (T2_fast_range(1)+T2_fast_range(2))/2 signalData(1)]; % SF_fast, T2_slow, T2_fast, SI_0
opts.TolFun = 1e-10;
opts.TolX = 1e-10;
opts.Upper = [1 T2_slow_range(2) T2_fast_range(2) signalData(1)*100]; % SF_fast, T2_slow, T2_fast, SI_0

% Fit model to data
[fitresult, gof] = fit([bData, teData], signalData, ft, opts);


end