%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dielectric Slab Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize MATLAB
close all; clc;
clear all; 

% Dimensions
% Slab is 12 inches thick surrounded by air on each side
d = 12 * 2.54; %cm thick
dc = d/100; %meters Our critical dimension in this case is the whole slab
rNz = ceil(d); %This Nz represents real world size
rNz = rNz + 2;  % We are going to add air on each side of the problem.

%Material Vectors Initialized at Air
rER = ones([1 rNz]);
rUR = ones([1 rNz]);

% Add our Slab materials to the model
rER(1:rNz) = 6;
rUR(1:rNz) = 2;

% Frequency

freq_start = 0; %DC
freq_end = 1e9; %1Ghz

NFREQ = freq_end / 10e6; %Frequencies every 100Mhz upto 10Ghz
FREQ = linspace(freq_start, freq_end, NFREQ); %FREQ List

FDTD1D( dc, dc, rER, rUR, -1, -1, FREQ, NFREQ, 1000, -1, 'HW#6-P2-Dielectric Slab' );











