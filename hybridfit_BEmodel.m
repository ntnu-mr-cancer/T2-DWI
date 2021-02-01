function [fitresult, gof] = hybridfit_BEmodel(b, signal)


% HYBRIDFIT_BEMODEL(B, SIGNAL)
% Creates a fit of signal fractions for DWI, using a bi-exponential model
% with a slow and a fast component.
%
% Model:
% SI/SI_0 = SF_slow exp(-b*ADC_slow) + SF_fast exp(-b*ADC_fast)
%         
% and:
% SF_slow + SF_fast = 1 
%
% Data for hybrid fit:
%	b: vector of b values [b1 b2]
%	signal: vector of pixel values [S_b1 S_b2]
% Output:
%	fitresult: a fit object representing the fit
%	gof: structure with goodness-of fit info


[bData, signalData] = prepareCurveData(b, signal);

% Specifying fixed ADC values for the components (values from literature)
% ADC_slow = 0.3; % µm²/ms
% ADC_fast = 2.6; % µm²/ms

% Setting up fittype and options
ft = fittype('SI_0*((1-SF_fast)*exp(-b*0.3)+SF_fast*exp(-b*2.6))','independent','b','dependent','SI','coefficients',{'SF_fast','SI_0'});
opts = fitoptions('Method','NonLinearLeastSquares');
opts.Display = 'Off';
opts.Lower = [0 0]; % SF_fast, SI_0
opts.MaxFunEvals = 10000;
opts.MaxIter = 1000;
opts.StartPoint = [0.5 signalData(1)]; % SF_fast, SI_0
opts.TolFun = 1e-10;
opts.TolX = 1e-10;
opts.Upper = [1 signalData(1)*100]; % SF_fast, SI_0

% Fit model to data
[fitresult, gof] = fit(bData, signalData, ft, opts);


end