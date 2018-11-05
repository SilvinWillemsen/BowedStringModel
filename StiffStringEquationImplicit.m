clear all;
close all;
clc;

fs = 44100;     % Sampling rate
gamma = 400;    % Wave-Speed
k = 1/fs;       % Time-step
s0 = 0.1;       % Damping coefficient
s1 = 0.005;

kappa = 2;      % Stiffness Factor

% Calculate grid spacing
h = sqrt((gamma^2*k^2 + 4 * s1 * k + sqrt((gamma^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
N = floor(1/h); % Number of gridpoints
h = 1/N; % Recalculate gridspacing

% Courant numbers
lambdaSq = (gamma*k/h)^2; 
muSq = (k * kappa / h^2)^2;

% Initialise state vectors
uPrev = zeros(N, 1);
u = zeros(N, 1);
uNext = zeros(N,1);

lengthSound = fs*5; % Set the length of the output
out = zeros(lengthSound, 1);

% Boundary condition (clamped, ss, free)
bc = "clamped";

%Bow Model
a = 100;                % free parameter
A = sqrt(2*a)*exp(1/2);

% User variables
Vb = 0.2;               % Bowing speed
Fb = 5000;               % Bowing force
bowPoint = floor(N/4);  % Bowing Position
pickup = floor(N/3);    % Pickup position

% Initialise variables for Newton Raphson 
tol = 1e-4;
qSave = zeros (lengthSound, 1);
qPrev2 = 0;

drawString = false;

for t = 1 : lengthSound
    if bc == "clamped"
        
        % Newton Raphson
        for l = 2 : N - 1
            if l > 2 && l < N - 1
                if l == bowPoint
                    b = 2/k * Vb - (2/k^2)*(u(l)-uPrev(l)) -(gamma/h)^2*(u(l+1) - 2*u(l) + u(l-1))...
                        + kappa^2/h^4 * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) + 2 * s0 * Vb ...
                        - (2*s1/(k*h^2)) * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                          - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)));
                    eps = 1;
                    i = 0;
                    while eps>tol
                        q=qPrev2-(Fb*A*qPrev2*exp(-a*qPrev2^2)+2*qPrev2/k+b)/...
                         (Fb*A*(1-2*a*qPrev2^2)*exp(-a*qPrev2^2)+2/k);
                        eps = abs(q-qPrev2);
                        qPrev2 = q;
                        i = i + 1;
                    end
                    if l == bowPoint
                        qSave(t) = q;
                    end
                    uNext(l) = 2*k*(q+Vb)+uPrev(l);
                else
                    uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                          - muSq * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) ...
                          + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                          * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                          - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                          / (1 + s0 * k);
                end
            end
        end
    end
    
    if drawString == true && mod(t,5) == 0
        plot(u, 'Linewidth', 1);
        drawnow;
    end
    out(t) = uNext(pickup);
    uPrev = u;
    u = uNext;
end

winSize = 2^12;
overlap = 2^8;
% spectrogram(out*100,hanning(winSize),overlap,[],fs, 'MinThresh', -100,'yaxis')
plot(out);