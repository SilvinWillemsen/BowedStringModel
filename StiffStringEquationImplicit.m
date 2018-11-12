clear all;
close all;
clc;

fs = 44100;     % Sampling rate
gamma = 400;    % Wave-Speed
k = 1/fs;       % Time-step
s0 = 0.1;       % Damping coefficients
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
uNext = zeros(N, 1);

uIPrev = zeros(N, 1);
uI = zeros(N, 1);
uINext = zeros(N, 1);

lengthSound = fs*5; % Set the length of the output
out = zeros(lengthSound, 1);

% Boundary condition (clamped, ss, free)
bc = "clamped";

%Bow Model
a = 100;                % free parameter
A = sqrt(2*a)*exp(1/2);

% User variables
Vb = 0.2;               % Bowing speed
Fb = 50;               % Bowing force
pickup = floor(N/3);    % Pickup position

% Initialise variables for Newton Raphson 
tol = 1e-4;
qSave = zeros (lengthSound, 1);
qPrev = 0;

drawString = false;

vectorRepresentation = true;

vec = 3:N-2;

% Introduce an extra gridpoint at the beginning of the interpolated state vectors 
uINext = zeros(N+1,1);
uI = zeros(N+1,1);
uIPrev = zeros(N+1,1);

% Set start and end time of the moving bow 
start = fs;
ending = fs*2;

% Set start and end position of the moving bow
startPoint = N/4 + 1;
endPoint = N/4 + 2;

for t = 1 : lengthSound
    if bc == "clamped"
        
        % Move the bow point
        if t < start
            bowPoint = startPoint;
        elseif t >= start && t < ending
            bowPoint = startPoint + (endPoint-startPoint) * (t-start) / (ending-start);
        else
            bowPoint = endPoint;
        end
        alpha = bowPoint - floor(bowPoint);

        % Linearly interpolate entire string
        uSave = u;
        uI(2:end-1) = u(1:end-1) * (1-alpha) + u(2:end) * alpha;
        uIPrev(2:end-1) = uPrev(1:end-1) * (1-alpha) + uPrev(2:end) * alpha;
        
        % Set the bowpoint
        bp = floor(bowPoint) + 1;
        
        J = 1/h;
        
        % Newton-Raphson
        b = 2/k * Vb - (2/k^2)*(uI(bp)-uIPrev(bp)) -(gamma/h)^2*(uI(bp+1) - 2*uI(bp) + uI(bp-1))...
            + kappa^2/h^4 * (uI(bp+2) - 4*uI(bp+1) + 6*uI(bp) - 4*uI(bp-1) + uI(bp-2)) + 2 * s0 * Vb ...
            - (2*s1/(k*h^2)) * ((uI(bp+1) - 2 * uI(bp) + uI(bp-1)) ...
              - (uIPrev(bp+1) - 2 * uIPrev(bp) + uIPrev(bp-1)));
        eps = 1;
        i = 0;
        qPrev=0;
        while eps>tol
            q=qPrev-(J*Fb*A*qPrev*exp(-a*qPrev^2)+2*qPrev/k+b)/...
             (J*Fb*A*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k);
            eps = abs(q-qPrev);
            qPrev = q;
            i = i + 1;
            if i > 10000
                disp("whut")
            end
        end
        qSave(t) = q;

        % Solve the FDE for the interpolated string
        vec = 3:length(uI)-2;
        uINext(vec) = (2 * uI(vec) - uIPrev(vec) + lambdaSq * (uI(vec+1) -2 * uI(vec) + uI(vec-1)) ...
              - muSq * (uI(vec+2) - 4*uI(vec+1) + 6*uI(vec) - 4*uI(vec-1) + uI(vec-2)) ...
              + s0 * k * uIPrev(vec) + (2*s1*k)/h^2 ...
              * ((uI(vec+1) - 2 * uI(vec) + uI(vec-1)) ...
              - (uIPrev(vec+1) - 2 * uIPrev(vec) + uIPrev(vec-1)))) ...
              / (1 + s0 * k);
        uINext(bp) = 2 * k * (q + Vb) + uIPrev(bp);

        % Interpolate back to the original grid point locations (and make sure that clamped boundary condition is honoured)
        uNext(3:N-2) = uINext(3:end-3) * alpha + uINext(4:end-2) * (1-alpha);
        u(3:N-2) = uI(3:end-3) * alpha + uI(4:end-2) * (1-alpha);
    end
    
    % Plot
    if drawString == true && mod(t,1) == 0 && t >= start - 100
        lowRange = 17;
        highRange = 27;
        range = lowRange:highRange;
        
        subplot(2,1,1);
        cla
        plot(u);
        title("Total string")
%         range = 1:N;
        yMin = -0.2e-4;
        yMax = 0.2e-4;
        ylim([yMin, yMax]);
        hold on;
        plot([lowRange-1 lowRange-1], [yMin yMax], 'k')
        plot([highRange highRange], [yMin yMax], 'k')
        subplot(2,1,2);
        cla
        plot(range, uSave(range), 'o-')
        hold on; 
        plot(range + alpha - 1, uI(range), '.-')
        plot(range, u(range), '+-')
        legend("Original", "Interpolated", "Result")
        title(strcat("\alpha = ", string(alpha)));
        drawnow;
    end
    
    out(t) = uNext(pickup);
    
    uPrev = u;
    u = uNext;
end

plot(out);