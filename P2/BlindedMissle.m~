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

% Frequency we want to transmit
lambda_0 = 980e-9 % nm;
f0 = c0/lamda_0

% Bragg Grating Materials
nSiN = 2.0;
erSiN = nSiN^2;

nSiO2 = 1.5;
erSiO2 = nSiO2^2;



%Calculate the Length of our layers.
LSiN = lambda_0/(4*nSiN);
LSiO2 = lambda_0/(4*nSiO2);

% We need to place a Anti-Reflective Layer on each side of the Radome
e_radome = 12;
e_air = 1;
e_nonreflective = sqrt(e_radome*e_air);

n_nonreflective = sqrt(e_nonreflective);
d_nonreflective = lambda_trans/(4*n_nonreflective);

dc = LSiN+LSiO2; %meters Our critical dimension in this case the anti-reflectivelayer
rNz = ceil((round(d_radome*100) + 2*round(d_nonreflective*100)))+2; %This Nz represents real world size

% Material Vectors Initialized at Air
rER = ones([1 rNz]);
rUR = ones([1 rNz]); 

% Add our Materials to the model
zstart = 1;
zend = ceil(d_nonreflective*100);
rER(zstart: zend) = e_nonreflective;

zstart = zend + 1;
zend = zstart + floor(d_radome*100);
rER(zstart:zend) = e_radome;

zstart=zend+1;
zend = zstart + floor(d_nonreflective*100);
rER(zstart:zend) = e_nonreflective;


% Frequency

freq_start = 0; %DC
freq_end = 5e9;%f_trans*2; %1Ghz

NFREQ = freq_end / 10e6; %Frequencies every 100Mhz upto 5Gz
FREQ = linspace(900, 1100, 100) * 1e-9; %FREQ List

FDTD1D( dc, (d_radome+2*d_nonreflective), rER, rUR, 35000, 100, FREQ, NFREQ, 2000);











