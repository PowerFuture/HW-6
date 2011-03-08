%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Invisible Slab Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
close all; clc;
clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Problem
%    A radome is being deigned to protect an antenna.
%    Antenna operates at 2.4Ghz
%    radome is 1ft thick with a dielectric constant = 12
%    We want  to maximize transmission through dome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
c0 = 299792458;

% Frequency we want to transmit
f_trans = 2.4e9; %2.4Ghz
lambda_trans = c0/f_trans; 


% Dimensions
% Radome is 12 inches thick
d_radome = 12 * 2.54/100; %cm thick

% We need to place a Anti-Reflective Layer on each side of the Radome
e_radome = 12;
e_air = 1;
e_nonreflective = sqrt(e_radome*e_air);

n_nonreflective = sqrt(e_nonreflective);
d_nonreflective = lambda_trans/(4*n_nonreflective);

dc = d_nonreflective; %meters Our critical dimension in this case the anti-reflectivelayer
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
FREQ = linspace(freq_start, freq_end, NFREQ); %FREQ List

FDTD1D( dc, (d_radome+2*d_nonreflective), rER, rUR, 35000, 100, FREQ, NFREQ, 2000);











