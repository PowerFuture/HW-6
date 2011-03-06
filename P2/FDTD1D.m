function FDTD1D( dc, rER, rUR, Steps, FREQ, NFREQ )
%FDTD1D Method executes a FDTD1D Model
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pre-Program Work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constants
c0 = 299792458; %m/s
e0 = 8.854187817*10^-12; %F/m
u0 = 1.256637061*10^-6; %H/m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization of Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_max = FREQ(length(FREQ));
nmax = Getnmax(rER, rUR);
%Compute Grid Resolution
% Wave Length Rsolution
N_lambda = GetNlambda(rER, rUR);
lambda_min = c0 / (f_max*nmax);
d_lambda = lambda_min/N_lambda;

% Structure Resolution
N_d = 4;
d_d = dc/4; 

% Calculate grid resolution dz
dz = min(d_lambda, d_d);
Nz = ceil(dc/dz);
dz = dc/Nz;

% Add free space buffer and TF/SF
buffer = ceil(d_lambda/dz) * 5;
buffert = buffer*2 + 3;
Nz = Nz + buffert;

%Compute Time Steps
dt = dz/(2*c0); %secs

% Source Parameters
nzc = 2;  %Position of Sources at our TF/SF boundary
tau = 0.5/f_max;        % tau parameter
t0 = 6*tau;              % Delay/Pulse Position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cf = floor((Nz - buffert)/length(rUR)); % Conversion factor to convert our real grid to our numerical grid

%Material Vectors
ER = zeros([1 Nz]);
UR = zeros([1 Nz]);
 
% We Need to lay our real materials vectors over our numerical material
% grid
 
% Lets place our real grid in proper location on numerical grid
for i = 0 : length(rER)-1
  index = buffer+2 + i*cf+1;
  %disp(['i: ' num2str(i) ' i2: ' num2str(index)]);
  ER(index) = rER(i+1);
  UR(index) = rUR(i+1);
end

% Need to backfill in our values
ER(1:buffer+2) = 1;
ER(length(ER)-buffer-1:length(UR)) = 1;
UR(1:buffer+2) = 1;
UR(length(UR)-buffer-1:length(UR)) = 1;
 
for i=buffer+2 : length(ER-buffer-1)
  if(ER(i) == 0)
    ER(i) = ER(i-1);
  end
  
  if(UR(i) == 0)
    UR(i) = UR(i-1);
  end
end

%ER(5) = 3.4641;
%ER(31) = 3.4641;
%ER(6:30) = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STEPS = Steps;
if(STEPS == -1)
  tprop = (nmax*Nz*dz)/c0; % Wave Propagation time;
  T = 12*tau + 5*tprop;
  STEPS = ceil(T/dt);
end

ta = [0:STEPS-1]*dt;     % Time Axis;

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

%Grid Axis
za=[0:Nz-1]*dz;


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


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Parameters');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

disp(['f_max' num2str(f_max)]);
disp(['lamda_min: ' num2str(lambda_min)]);
disp(['d_lambda: ' num2str(d_lambda)]);
disp(['nmax: ' num2str(nmax)]);
disp(['dc: ' num2str(dc)]);
disp(['d_d: ' num2str(d_d)]);
disp(['Nz: ' num2str(Nz)]);
disp(['buffer: ' num2str(buffer)]);
disp(['dz: ' num2str(dz)]);
disp(['dt: ' num2str(dt)]);
disp(['tau: ' num2str(tau)]);
disp(['t0: ' num2str(t0)]);
disp(['STEPS: ' num2str(STEPS)]);
disp(['s: ' num2str(s)]);
disp(['A: ' num2str(A)]);
disp(ER);

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
 
 
 if(mod(t,20) == 0)
   h = subplot(11,1,1:4);
   Draw1D(ER, Ey, Hx, dz);
   axis([za(1) za(Nz) -1.5 1.5]);
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








end

