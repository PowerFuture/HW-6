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

nanometers = 1e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulated Environment Settings
NLAM = 5e9;  % 5Ghz
lambda_0 = 980*nanometers; % Anti-reflective frequency;
lambda_min = 900*nanometers;
PERIODS = 15;

nSiN = 2.0;
erSiN = nSiN^2;

nSiO2 = 1.5;
erSiO2 = nSiO2^2;
nmax = nSiN;

%Calculate the Length of our layers.
LSiN = lambda_0/(4*nSiN); 
LSiO2 = lambda_0/(4*nSiO2);

dc = LSiN;%meters Our critical dimension in this case is the width of one period.


%Compute Grid Resolution
N_lambda = 20;
d_wl = lambda_min/N_lambda/nmax;
N_d = 4;
d_d = dc/4; % since we are only working with freespace we will set d to 1;
dz = min(d_wl, d_d);
Nprime = ceil(dc/dz);
dz = dc/Nprime;

Nz = PERIODS*ceil((LSiN+LSiO2)/dz);


Nz = Nz + 2*(100) + 3;


%Grid Axis
za=[0:Nz-1]*dz;

%Compute Time Steps
dt = dz/(2*c0); %secs

% Source Parameters
nzc = 2;  %Position of Sources
NLAM = 100;
LAMBDA = linspace(900, 1100, NLAM)*nanometers; %FREQ List
tau = 0.5/(c0/lambda_min);        % tau parameter
t0 = 6*tau;              % Delay/Pulse Position

T = 12*tau + 5*(nmax*Nz*dz/c0);
STEPS = ceil(T/dt);
STEPS=20000;

ta = [0:STEPS-1]*dt;     % Time Axis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material Vectors
ER = ones([1 Nz]);
UR = ones([1 Nz]);

nstart = 100+2;
nend = nstart;


for p = 1:PERIODS
  nend = nstart + round(LSiN/dz)-1;
  ER(nstart:nend) = erSiN;
  
  nstart = nend+1;
  nend = nstart+ round(LSiO2/dz) -1;
  ER(nstart:nend) = erSiO2;
   
  nstart = nend+1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = dz/(2*c0) + dt/2;    % Delay between E and H
Esrc = exp(-((ta-t0)/tau).^2); % E Source
A = -sqrt(ER(nzc)/UR(nzc));    % H Amplitude
Hsrc = A*exp(-((ta-t0+s)/tau).^2); % H Source



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
REF = zeros(1, NLAM);
TRN = zeros(1, NLAM);
SRC = zeros(1, NLAM);
K = exp(-1i*2*pi*dt*(c0./LAMBDA));

SSFK = exp(-1i*2*pi*dt*c0./lambda_0);
SSFPOWER = zeros(1, Nz);
SSFSRC = zeros(1, Nz);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Parameters');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

disp(['lamda_min: ' num2str(lambda_min)]);
disp(['d_lambda: ' num2str(d_wl)]);
disp(['nmax: ' num2str(nmax)]);
disp(['dc: ' num2str(dc)]);
disp(['d_d: ' num2str(d_d)]);
disp(['Nz: ' num2str(Nz)]);
disp(['dz: ' num2str(dz)]);
disp(['Length: ' num2str(Nz*dz)]);
disp(['dt: ' num2str(dt)]);
disp(['tau: ' num2str(tau)]);
disp(['t0: ' num2str(t0)]);
disp(['STEPS: ' num2str(STEPS)]);
disp(['s: ' num2str(s)]);
disp(['A: ' num2str(A)]);
% disp(['ER: ' num2str(length(ER))]);
% disp(ER);
% disp(['UR: ' num2str(length(UR))]);;
% disp(UR);
% return;



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
 for nf = 1: NLAM
   REF(nf) = REF(nf) + (K(nf)^t)*Ey(1)*dt;
   TRN(nf) = TRN(nf) + (K(nf)^t)*Ey(Nz)*dt;
   SRC(nf) = SRC(nf) + (K(nf)^t)*Esrc(t)*dt;
 end
 
  for n = 3 : Nz-1
   SSFPOWER(n) = SSFPOWER(n) + (SSFK^t)*Ey(n)*dt;
   SSFSRC(n) = SSFSRC(n) + (SSFK^t)*Esrc(t)*dt;
  end
   
 if(mod(t,1000) == 0)
   h = subplot(11,1,1:4);
   Draw1D(ER, Ey, Hx, dz);
   axis([za(1) za(Nz) -1.5 1.5]);
   xlabel('z');
   title(['Field at Step ' num2str(t) ' of ' num2str(STEPS)]);    
 
   R = abs(REF./SRC).^2;
   T = abs(TRN./SRC).^2;

   subplot(11,1,8:11)
   plot(LAMBDA/nanometers, 10*log10(R), '-r'); hold on;
   plot(LAMBDA/nanometers, 10*log10(T), '-b');
   plot(LAMBDA/nanometers, 10*log10(R+T), ':k', 'LineWidth', 2); hold off;
   axis([LAMBDA(1)/nanometers LAMBDA(NLAM)/nanometers -50 0]);
   xlabel('Frequency');
   title('Reflectance and Transmittance');
   drawnow();
 end
   

  
            
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

SSFPOWER = abs(SSFPOWER./SSFSRC).^2;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
SetFigure(fig, 'HW#6-P2-Blinded Missle', [500 274 965 826]);

subplot(11,1,1:4);
plot(LAMBDA/nanometers, 10*log10(R), '-r'); hold on;
plot(LAMBDA/nanometers, 10*log10(T), '-b', 'LineWidth', 2);
plot(LAMBDA/nanometers, 10*log10(R+T), ':k', 'LineWidth', 2); hold off;
axis([LAMBDA(1)/nanometers LAMBDA(NLAM)/nanometers -50 0]);
xlabel('Frequency');
title('Reflectance and Transmittance');
legend('Reflectance','Transmittance','Consersvation','Location','SouthEast');


subplot(11,1,8:11)
plot(SSFPOWER, '-r', 'LineWidth', 2);
axis([3 Nz-1  -0.1  5]);
xlabel('Nz');
title('Steady State Field');
