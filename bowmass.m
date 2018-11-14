% matlab script bowmass.m
% finite difference scheme for a bowed mass-spring system
% soft friction characteristic w/iterative Newton-Raphson method
%%%%%% begin global parameters
SR = 44100; % sample rate (Hz)
f0 = 200; % oscillator frequency (Hz)
FB = 500; % bow force/mass (m/s^2)
TF = 0.1; % simulation duration (s)
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
f = zeros(NF,1);
vr = zeros(NF,1);
qlast = 0;
% time step restrictions
if(k>min(1/(pi*f0),exp(1)/(FB*sqrt(2*sig))))
    error("Time step too large");
end

%%%%%% start main loop

for n=3:NF
    % Newton-Raphson method to determine relative velocity
    b = (2*pi*f0)^2*u(n-1)-(2/k^2)*(u(n-1)-u(n-2))+(2/k)*vB;
    eps = 1;
    while eps>tol
        q=qlast-(FB*A*qlast*exp(-sig*qlast^2)+2*qlast/k+b)/...
         (FB*A*(1-2*sig*qlast^2)*exp(-sig*qlast^2)+2/k);
        eps = abs(q-qlast);
        qlast = q;
    end
    % update position of mass and relative bow velocity
    u(n) = 2*k*(q+vB)+u(n-2);
    vr(n-1) = q;
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