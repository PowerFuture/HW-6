%FDTD1D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pre-Program Work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
close all; clc;
clear all; 

%Constants
c0 = 299792458; %m/s
e0 = 8.854187817*10^-12; %F/m
u0 = 1.256637061*10^-6; %H/m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulated Environment Settings
STEPS = 2000;
Nz = 180;
dz = 0.02;
f_max = 1e9;  % 10Ghz
dc = .3048;
nmax = 3.46;

%Compute Grid Resolution
N_lambda = 20;
lambda_min = c0 / (f_max*nmax);
d_wl = lambda_min/N_lambda;
N_d = 4;
d_d = dc/4; % since we are only working with freespace we will set d to 1;
dz = min(d_wl, d_d);
Nz = ceil(dc/dz);
dz = dc/Nz;

Nz = Nz + 2*(10) + 3;

disp(['lambda_min: ' num2str(lambda_min)]);
disp(['d_wl: ' num2str(d_wl)]);
disp(['dc: ' num2str(dc)]);
disp(['d_d: ' num2str(d_d)]);
disp(['Nz: ' num2str(Nz)]);
disp(['dz: ' num2str(dz)]);

%Grid Axis
za=[0:Nz-1]*dz;

%Compute Time Steps
dt = dz/(2*c0); %secs

% Source Parameters
nzc = round (2);  %Position of Sources
NFREQ = f_max / 10e6; %Frequencies every 100Mhz upto 10Ghz
FREQ = linspace(0, f_max, NFREQ); %FREQ List
tau = 0.5/f_max;        % tau parameter
t0 = 6*tau;              % Delay/Pulse Position

T = 12*tau + 5*(nmax*Nz*dz/c0);
STEPS = ceil(T/dt);


disp(['dt: ' num2str(dt)]);
disp(['tau: ' num2str(tau)]);
disp(['t0: ' num2str(t0)]);
disp(['Steps: ' num2str(STEPS)]);


ta = [0:STEPS-1]*dt;     % Time Axis;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material Vectors
ER = ones([1 Nz]);
UR = ones([1 Nz]);

nz1 = 2+10+1;
nz2 = nz1 + round(dc/dz) - 1;

disp(['nz1: ' num2str(nz1)]);
disp(['nz2: ' num2str(nz2)]);

UR(nz1:nz2) = 2;
ER(nz1:nz2) = 6;

disp(ER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = dz/(2*c0) + dt/2;    % Delay between E and H
Esrc = exp(-((ta-t0)/tau).^2); % E Source
A = -sqrt(ER(nzc)/UR(nzc));    % H Amplitude
Hsrc = A*exp(-((ta-t0+s)/tau).^2); % H Source

disp(['s: ' num2str(s)]);
disp(['A: ' num2str(A)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FDTD Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Update Coefficients
mER = (c0*dt/dz)./ER;
mHR = (c0*dt/dz)./UR;

% Initialize Feilds
Ey = zeros([1 Nz]);
Hx = zeros([1 Nz]);


%PAB Parameters
h1 = 0; h2 = 0; h3 = 0;
e1 = 0; e2 = 0; e3 = 0;

%Power Measurements
REF = zeros(1, NFREQ);
TRN = zeros(1, NFREQ);
SRC = zeros(1, NFREQ);
K = exp(-1i*2*pi*dt*FREQ);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:STEPS
 
  % Calculate H
  for nz = 1:Nz-1
    Hx(nz) = Hx(nz) + mHR(nz)*(Ey(nz+1)-Ey(nz));
  end
  
  Hx(Nz) = Hx(Nz) + mHR(Nz)*(e3 - Ey(Nz));

  %H Sources
  Hx(nzc-1) = Hx(nzc-1) - mHR(nzc-1)*Esrc(t);

  h3 = h2; h2 = h1; h1 = Hx(1); % Boundary Params;
  
  % Calculate E  
  Ey(1) = Ey(1) + mER(1)*(Hx(1) - h3);
  for nz = 2:Nz
    Ey(nz) = Ey(nz) + mER(nz)*(Hx(nz)-Hx(nz-1)); 
  end
  
  %Inject Source
  Ey(nzc) = Ey(nzc) - mER(nzc)*Hsrc(t);

  e3=e2; e2=e1; e1=Ey(Nz); % Boundary Params;
 
 %Update Fourier Transforms
 for nf = 1: NFREQ
   REF(nf) = REF(nf) + (K(nf)^t)*Ey(1)*dt;
   TRN(nf) = TRN(nf) + (K(nf)^t)*Ey(Nz)*dt;
   SRC(nf) = SRC(nf) + (K(nf)^t)*Esrc(t)*dt;
 end
 
 
 if(mod(t,10) == 0)
   h = subplot(11,1,1:4);
   Draw1D(ER, Ey, Hx, dz);
   axis([za(1) za(Nz) -1.1 1.1]);
   xlabel('z');
   title(['Field at Step ' num2str(t) ' of ' num2str(STEPS)]);    
 
   R = abs(REF./SRC).^2;
   T = abs(TRN./SRC).^2;

   subplot(11,1,8:11)
   plot(FREQ, R, '-r'); hold on;
   plot(FREQ, T, '-b');
   plot(FREQ, R+T, ':k', 'LineWidth', 2); hold off;
   axis([FREQ(1) FREQ(NFREQ) -0.1 1.5]);
   xlabel('Frequency');
   title('Reflectance and Transmittance');
 end
   
drawnow();
  
            
  %if(mod(t,50) == 0)
  %  saveas(h, ['images/' num2str(t) '.jpg'], 'jpg');
  %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REF = abs(REF./SRC).^2;
TRN = abs(TRN./SRC).^2;
CON = REF+TRN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
SetFigure(fig, 'HW#5-P2', [500 274 965 826]);

plot(FREQ, REF, '-r', 'LineWidth', 2); hold on;
plot(FREQ, TRN, '-b', 'LineWidth', 2);
plot(FREQ, CON, ':k', 'LineWidth', 3); hold off;
axis([FREQ(1) FREQ(NFREQ) -0.1 1.2]);
xlabel('Frequency');
title('Reflectance and Transmittance');




