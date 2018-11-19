% matlab script bowmass.m
% finite difference scheme for a bowed mass-spring system
% soft friction characteristic w/iterative Newton-Raphson method

clear all;
close all;
clc;

%%%%%% begin global parameters
SR = 44100; % sample rate (Hz)
f0 = 200; % oscillator frequency (Hz)
FB = 500; % bow force/mass (m/s^2)
TF = 1; % simulation duration (s)
vB = 0.2; % bow velocity (m/s)
sig = 100; % friction law free parameter (1/m^2) 
tol = 1e-4; % tolerance for Newton-Raphson method
%%%%%% end global parameters


% derived parameters
NF = floor(TF*SR);
k = 1/SR;
A = exp(1/2)*sqrt(2*sig);

% initialize time series/iterative method

u = zeros(NF,1);
u2 = zeros(NF,1);
f = zeros(NF,1);
vr = zeros(NF,1);
qlast = 0;
qlast2 = 0;
% time step restrictions
if(k>min(1/(pi*f0),exp(1)/(FB*sqrt(2*sig))))
    error("Time step too large");
end
%%%%%% start main loop

for n=3:NF
    % Newton-Raphson method to determine relative velocity
    b = (2*pi*f0)^2*u(n-1)-(2/k^2)*(u(n-1)-u(n-2))+(2/k)*vB;
    eps = 1;
    i = 0;
    while eps>tol
        q=qlast-(FB*A*qlast*exp(-sig*qlast^2)+2*qlast/k+b)/...
         (FB*A*(1-2*sig*qlast^2)*exp(-sig*qlast^2)+2/k);
        eps = abs(q-qlast);
        qlast = q;
        i = i + 1;
        if i > 10000
            disp('whut')
        end
    end
    
    b2 = (2*pi*f0)^2*u2(n-1)-(2/k^2)*(u2(n-1)-u2(n-2))+(2/k)*vB;
    eps = 1;
    i = 0;
    while eps>tol
        q2=qlast2-(FB*A*qlast2*exp(-sig*qlast2^2)+2*qlast2/k+b2)/...
         (FB*A*(1-2*sig*qlast2^2)*exp(-sig*qlast2^2)+2/k);
        eps = abs(q2-qlast2);
        qlast2 = q2;
        i = i + 1;
        if i > 10000
            disp('whut')
        end
    end
    % update position of mass and relative bow velocity
    u(n) = 2*u(n-1) - u(n-2) - k^2*(2*pi*f0)^2 * u(n-1) - k^2 * FB * sqrt(2*sig) * q * exp(-sig*q^2 + 1/2);
    u2(n) = 2*k*(q2+vB)+u2(n-2);
%     vr(n-1) = q;
end
%%%%%% end main loop

% plot mass displacement and relative bow velocity
tax = [0:NF-1]*k;
subplot(2,1,1); plot(tax, u, "k"); title("Displacement of Mass");
xlabel("time (s)"); 
subplot(2,1,2); 
plot(tax, vr, "k");
title("Relative Bow Velocity");
xlabel("time (s)");