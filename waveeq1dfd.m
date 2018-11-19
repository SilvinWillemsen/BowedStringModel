% matlab script waveeq1dfd.m
% finite difference scheme for the 1D wave equation % fixed boundary conditions
% raised cosine initial conditions

%%%%%% begin global parameters

clear all; 
close all;

SR = 44100;             % sample rate (Hz)
f0 = 125;               % fundamental frequency (Hz)
TF = 5;                 % duration of simulation (s)
ctr = 0.7; wid = 0.1;   % center location/width of excitation
u0 = 1; v0 = 0;         % maximum initial displacement/velocity 
rp = 0.3;               % position of readout (0-1)
lambda = 1;             % Courant number

%%%%%% end global parameters

% begin derived parameters
gamma = 2*f0;       % wave equation free parameter
k = 1/SR;           % time step
NF = floor(SR*TF);  % duration of simulation (samples)

% stability condition/scheme parameters
h = gamma*k/lambda; 
N = floor(1/h); 
h = 1/N; 
lambda = gamma*k/h;
s0 = 2*(1-lambda^2); s1 = lambda^2;

% readout interpolation parameters
rp_int = 1+floor(N*rp);     % rounded grid index for readout
rp_frac = 1+rp/h-rp_int;    % fractional part of readout location

% create raised cosine
xax = [0:N]'*h;
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));

% initialize grid functions and output

% u2 = u0*rc;
% u1 = (u0+k*v0)*rc;
u2 = zeros(N+1,1);
u1 = zeros(N+1,1);
u = zeros(N+1,1);
uPrev = zeros(N+1,1);
uCur = zeros(N+1,1);
uNext = zeros(N+1,1);
out = zeros(NF,1); 

bowPoint = floor(0.3*N);
bowVector = zeros(length(u), 1);
bowVector(bowPoint) = 1;

Fb = 10;
Vb = 0.2;
a = 100;
A = sqrt(2*a)*exp(1/2);
%%%%%% start main loop
tol = 1e-4;
qPrev = 0;
qPrev2 = 0;
vectMode = false;
i = 0;
for n=3:NF
    % Vector representation
    I = zeros(N,1);
    I(bowPoint) = 1;
        b = -(gamma/h)^2*(u1(bowPoint + 1) - 2*u1(bowPoint) + u1(bowPoint-1)) - (2/k^2)*(u1(bowPoint)-u2(bowPoint))+(2/k)*Vb;
        eps = 1;
        qPrev = 0;
        while eps>tol
            q=qPrev-(1/h*Fb*A*qPrev*exp(-a*qPrev^2)+2*qPrev/k+b)/...
             (1/h*Fb*A*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k);
            eps = abs(q-qPrev);
            qPrev = q;
            i = i + 1;
            if i > 10000000
                disp("whut")
            end
        end 
        iSave(n) = i;
        qSave(n) = q;
%         excitation = 1/h*I(2:N)*Fb*A*exp(-a*q^2);
        u(2:N) = -u2(2:N)+s0*u1(2:N)+s1*(u1(1:N-1)+u1(3:N+1));%-excitation; % scheme calculation
        u(bowPoint) = 2*k*(q+Vb)+u2(bowPoint);
        out(n) = (1-rp_frac)*u(rp_int)+rp_frac*u(rp_int+1);
        u2 = u1; u1 = u; % update of grid variables
        
%         %For-loop Representation
%         for l = 2:N
%             if l == bowPoint
%                 b = -(gamma/h)^2*(uCur(l+1) - 2*uCur(l) + uCur(l-1)) - (2/k^2)*(uCur(l)-uPrev(l))+(2/k)*Vb;
%                 eps = 1;
%                 while eps>tol
%                     q=qPrev2-(Fb*A*qPrev2*exp(-a*qPrev2^2)+2*qPrev2/k+b)/...
%                      (Fb*A*(1-2*a*qPrev2^2)*exp(-a*qPrev2^2)+2/k);
%                     eps = abs(q-qPrev2);
%                     qPrev2 = q;
%                 end 
%                 uNext(l) = 2*k*(q+Vb)+uPrev(l);
%             else
%                 uNext(l) = s0 * uCur(l)+ s1 * (uCur(l+1) + uCur(l-1)) - uPrev(l);
%             end
%         end
%         out2(n) = (1-rp_frac)*uNext(rp_int)+rp_frac*uNext(rp_int+1);
%         uPrev = uCur;
%         uCur = uNext;
%         hold on;
%         plot(uNext); drawnow;
%     end
end
% readout
%%%%%% end main loop
% plot output waveform
plot([0:NF-1]*k, out, 'k');
xlabel('t'); ylabel('u'); title('1D Wave Equation: FD Output'); axis tight
hold on; plot(out2);
% play sound
% soundsc(out,SR);