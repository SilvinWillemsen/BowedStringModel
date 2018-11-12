clear all;
close all;
clc;

fs = 44100;     % Sampling rate
f0 = 196.00;    % G3
gamma = f0*2;     % Wave-Speed
k = 1/fs;       % Time-step
s0 = 0.1;       % Damping coefficients
s1 = 0.005;

B = 0.001; %inharmonicity coefficient
kappa = sqrt(B)*(gamma/pi); % Stiffness Factor     

% Calculate grid spacing
h = sqrt((gamma^2*k^2 + 4 * s1 * k + sqrt((gamma^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
N = floor(1/h); % Number of gridpoints
h = 1/N; % Recalculate gridspacing

% Courant numbers
lambdaSq = (gamma*k/h)^2; 
muSq = (k * kappa / h^2)^2;

% Initialise state vectors
uPrev = zeros(N,1);
u = zeros(N,1);
uNext = zeros(N,1);  

lengthSound = fs*5; % Set the length of the output
out = zeros(lengthSound, 1);

% Boundary condition (clamped, ss, free)
bc = "clamped";

%Bow Model
a = 100;                % free parameter
BM = sqrt(2*a)*exp(1/2);

% User variables
Vb = 0.2;               % Bowing speed
Fb = 50;
pickup = floor(N/3);    % Pickup position

% Initialise variables for Newton Raphson 
tol = 1e-4;
qSave = zeros (lengthSound, 1);
qPrev = 0;

drawString = false;

% Matrices FDE
phi = (2 * s1 * k) / (h^2);
A = sparse(1:N, 1:N, (1 + s0*k) * ones(1, N), N, N);
B = (sparse(3:N, 1:N-2, -muSq * ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, (4 * muSq + phi + lambdaSq) * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (-6*muSq - 2 * phi - 2 * lambdaSq + 2) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, (4 * muSq + phi + lambdaSq) * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, -muSq * ones(1, N-2), N, N));
C = sparse(2:N, 1:N-1, -phi * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (s0 * k + 2 * phi - 1) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, -phi * ones(1, N-1), N, N);
            
% Matrices for b in NR 
kOh = (kappa/h)^2;
gOh = (gamma/h)^2;
phi = (2 * s1)/(k*h^2);
bB = (sparse(3:N, 1:N-2, kOh * ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (6*kOh + 2 * gOh - 2/k^2 + 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, kOh * ones(1, N-2), N, N));
bC = sparse(2:N, 1:N-1, phi * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (2/k^2 - 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, phi * ones(1, N-1), N, N);

% Set start and end time of the moving bow 
start = fs;
ending = fs*2;

% Set start and end position of the moving bow
startPoint = N/4;
endPoint = 3*N/4;

interpolation = "none"; % (none, linear, cubic) 
useExcitation = false; % Has to be true when using linear or cubic interpolation

for t = 1 : lengthSound
    
    % Move the bow point
    if t < start
        bowPoint = startPoint;
    elseif t >= start && t < ending
        bowPoint = startPoint + (endPoint-startPoint) * (t-start) / (ending-start);
    else
        bowPoint = endPoint;
    end
    
    % Set alpha (0-1)
    alpha = bowPoint - floor(bowPoint);

    % Set the bowpoint
    bp = floor(bowPoint);
    I = zeros(N,1);
    
    if interpolation == "linear"
        % linear interpolation
        I(bp) = (1-alpha);
        I(bp + 1) = alpha;
    elseif interpolation == "cubic"
        % cubic interpolation
        I(bp - 1) = (alpha * (alpha - 1) * (alpha - 2)) / -6;
        I(bp) = ((alpha - 1) * (alpha + 1) * (alpha - 2)) / 2;
        I(bp + 1) = (alpha * (alpha + 1) * (alpha - 2)) / -2;
        I(bp + 2) = (alpha * (alpha + 1) * (alpha - 1)) / 6;
    elseif interpolation == "none"
        I(bp) = 1;
    end
    
    J = 1/h * I;

    % Newton-Raphson
    b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
    eps = 1;
    i = 0;
    while eps>tol
        q=qPrev-(1/h*Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0*qPrev+b)/...
         (1/h*Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0);
        eps = abs(q-qPrev);
        qPrev = q;
        i = i + 1;
        if i > 10000
            disp('whut');
        end
    end
    qSave(t) = q;

    if useExcitation == true 
        excitation = h^5*J*Fb*BM*exp(-a*q^2); % Need to scale with h^5 to get close to a bowing sound
        uNext = A\(B * u + C * uPrev) - excitation;
    else
    % Update uNext at the bowpoint as is done in A.3 of book
        uNext = A\(B * u + C * uPrev);
        uNext(bp) = 2 * k * (q + Vb) + uPrev(bp);
    end
    

    % Plot
    if drawString == true && mod(t,1) == 0 %&& t >= start - 100
        clf;
        plot(u);
        hold on;
        drawnow;
    end
    
    % Save output
    out(t) = uNext(pickup);
    
    % update state vectors
    uPrev = u;
    u = uNext;
end

plot(out);