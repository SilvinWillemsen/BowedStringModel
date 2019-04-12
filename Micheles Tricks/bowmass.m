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
TF = 0.1; % simulation duration (s)
vB = 0.1; % bow velocity (m/s)
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
u3 = zeros(NF,1);
f = zeros(NF,1);
vr = zeros(NF,1);
qlast = 0;
qlast2 = 0;
% time step restrictions
if(k>min(1/(pi*f0),exp(1)/(FB*sqrt(2*sig))))
    error("Time step too large");
end
%%%%%% start main loop
q2Prev = -vB;
q3Prev = -vB;
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
    qSave(n) = q;
    u(n) = 2*u(n-1) - u(n-2) - k^2*(2*pi*f0)^2 * u(n-1) - k^2 * FB * sqrt(2*sig) * q * exp(-sig*q^2 + 1/2);
    
%     q2 = (u2(n-1) - u2(n-2)) / k - vB;
    q2 = 2 / k * (u2(n-1) - u2(n-2)) - q2Prev - 2 * vB;
%     q3 = 2 / k * (u3(n-1) - u3(n-2)) - q3Prev - 2 * vB;
    q3 = (u3(n-1) - u3(n-2)) / k - vB;
%     q3Save(n) = q3;
    % update position of mass and relative bow velocity
    
    testFunc = FB * sqrt(2*sig) * exp(-sig * q2^2+1/2);
    u2(n) = (2 * u2(n-1) - u2(n-2) -k^2*(2*pi*f0)^2 * u2(n-1)  + testFunc * k / 2 * u2(n-2) + testFunc * k^2 * vB) / (1 + testFunc * k / 2);
    
    
    testFunc2 = k^2 * FB * sqrt(2*sig) * q3 * exp(-sig * q3^2+1/2);
    u3(n) = 2 * u3(n-1) - u3(n-2) - k^2*(2*pi*f0)^2 * u3(n-1) - testFunc2;
    
    q2Prev = q2;
    q3Prev = q3;
    q2Save(n) = q2;
    q3Save(n) = q3;
end
%%%%%% end main loop

% plot mass displacement and relative bow velocity
% tax = [0:NF-1]*k;
% subplot(2,1,1); plot(tax, u, "k"); title("Displacement of Mass");
% xlabel("time (s)"); 
% subplot(2,1,2); 
% plot(tax, vr, "k");
% title("Relative Bow Velocity");
% xlabel("time (s)");
timeVecUnderSamp = (0:length(u) - 1) / SR;
% timeVec = (0:length(u) - 1) / SR;
subplot(3,1,1)
plot(timeVecUnderSamp, u); hold on; plot(timeVecUnderSamp, u2); plot(timeVecUnderSamp, u3)
load saveTest.mat
plot(timeVec, saveNR)
legend(["NR","TRICKS","Backwards", "OversampledNR"])
subplot(3,1,[2,3])
range = 2600:2750;
plot(timeVecUnderSamp(range), u(range)); 
hold on; 
plot(timeVecUnderSamp(range), u2(range)); 
plot(timeVecUnderSamp(range), u3(range));
plot(timeVec(range*10), saveNR(range*10))
freqAxOS = ((0:length(saveNR) - 1) * SR * 10) / length(saveNR);

freqAx = ((0:length(u) - 1) * SR) / length(u);

legend(["NR","TRICKS","Backwards", "OversampledNR"])