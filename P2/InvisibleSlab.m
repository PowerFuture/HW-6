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

% Frequency
f_trans = 2.4e9; %2.4Ghz
lambda_trans = f_trans/c0;
% Dimensions
% Radome is 12 inches thick
d_radome = 12 * 2.54; %cm thick

% We need to place a Anti-Reflective Layer on each side of the Radome
e_radome = 12;
e_air = 1;
e_nonreflective = sqrt(e_radome*e_air);

n_nonreflective = sqrt(e_nonreflective);
d_nonreflective = lambda_trans/(4*n_nonreflective);

dc = d_nonreflective/100; %meters Our critical dimension in this case the anti-reflectivelayer
rNz = ceil(d_radome + 2*d_nonreflective); %This Nz represents real world size
rNz = rNz + 2;  % We are going to add air on each side of the problem.

% Material Vectors Initialized at Air
rER = ones([1 rNz]);
rUR = ones([1 rNz]);  % I don't think we change this value at all....

% Add our Materials to the model
zstart = 1;
zend = floor(d_nonreflective);
rER(zstart: zend) = e_nonreflective;

zstart = zend + 1;
zend = zstart + ceil(d_radome)+1;
rER(zstart:zend) = e_radome;

zstart=zend+1;
zend = zstart;
rER(zstart:zend) = e_nonreflective;



% Frequency

freq_start = 0; %DC
freq_end = 5e9;%f_trans*2; %1Ghz

NFREQ = freq_end / 10e6; %Frequencies every 100Mhz upto 5Gz
FREQ = linspace(freq_start, freq_end, NFREQ); %FREQ List

FDTD1D( dc, rER, rUR, -1, FREQ, NFREQ );











