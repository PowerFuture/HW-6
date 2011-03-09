%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Invisible Slab Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
close all; clc;
clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Problem
%    A missle is vulnerable to jamming from high power lasers
%    at lambda0 = 980nm.  Design a multilayer cover that would
%    prevent this energy from reachign the infrared camera.  
%    We need to provid at least 30db of suppression at 980 nm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
c0 = 299792458;

%units
nanometers = 1e-9;

% Frequency we want to transmit
lambda_0 = 980 * nanometers; % nm;
f0 = c0/lambda_0;

% Bragg Grating Materials
nSiN = 2.0;
erSiN = nSiN^2;

nSiO2 = 1.5;
erSiO2 = nSiO2^2;

Periods = 10;

%Calculate the Length of our layers.
LSiN = lambda_0/(4*nSiN) / nanometers; 
LSiO2 = lambda_0/(4*nSiO2) / nanometers;


dc = LSiN+LSiO2;%meters Our critical dimension in this case is the width of one period.
Length = ((LSiN+LSiO2)) * Periods+2;

rNz = ceil(Length); %This Nz represents real world size

% Material Vectors Initialized at Air
rER = ones([1 rNz]);
rUR = ones([1 rNz]); 

% Add our Materials to the model
nstart = 100+2;
nend = nstart;

disp(['Periods: ' PERIODS]);

for n = 1: PERIODS
  disp(['Period: ' num2str(n)]);
  nend = nstart + round(LSiN)-1;
  rER(nstart:nend) = erSiN;
  
  nstart = nend+1;
  nend = nstart+ round(LSiO2) -1;
  rER(nstart:nend) = erSiO2;
  
  nstart = nend+1;
end

% Frequency

freq_start = 900*nanometers/c0; %DC
freq_end = 1100*nanometers/c0;%f_trans*2; %1Ghz

NFREQ = freq_end / 100/c0; %Frequencies every 100nm
FREQ = linspace(freq_start, freq_end, NFREQ); %FREQ List

FDTD1D( dc*nanometers, Length*nanometers, rER, rUR, -1, -1, FREQ, NFREQ, 2000);











