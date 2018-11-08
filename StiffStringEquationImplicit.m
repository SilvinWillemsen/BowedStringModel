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
bowPoint = N/4;  % Bowing Position
pickup = floor(N/3);    % Pickup position

% Initialise variables for Newton Raphson 
tol = 1e-4;
qSave = zeros (lengthSound, 1);
qPrev = 0;
qPrev2 = 0;

drawString = true;

vectorRepresentation = true;

vec = 3:N-2;

uINext = zeros(N+1,1);
uI = zeros(N+1,1);
uIPrev = zeros(N+1,1);

% Set start/end time/position of the moving bow 
start = fs;
startPoint = N/4 + 1;
ending = fs*2;
endPoint = N/4 + 2;

for t = 1 : lengthSound
    if bc == "clamped"
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
        if alpha == 0
            bp = bowPoint + 1;
        else
            bp = ceil(bowPoint);
        end
        
        % Newton-Raphson
        b = 2/k * Vb - (2/k^2)*(uI(bp)-uIPrev(bp)) -(gamma/h)^2*(uI(bp+1) - 2*uI(bp) + uI(bp-1))...
            + kappa^2/h^4 * (uI(bp+2) - 4*uI(bp+1) + 6*uI(bp) - 4*uI(bp-1) + uI(bp-2)) + 2 * s0 * Vb ...
            - (2*s1/(k*h^2)) * ((uI(bp+1) - 2 * uI(bp) + uI(bp-1)) ...
              - (uIPrev(bp+1) - 2 * uIPrev(bp) + uIPrev(bp-1)));
        eps = 1;
        i = 0;
        qPrev2=0;
        while eps>tol
            q=qPrev2-(1/h*Fb*A*qPrev2*exp(-a*qPrev2^2)+2*qPrev2/k+b)/...
             (1/h*Fb*A*(1-2*a*qPrev2^2)*exp(-a*qPrev2^2)+2/k);
            eps = abs(q-qPrev2);
            qPrev2 = q;
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

        % Interpolate back (and make sure that clamped boundary condition is honoured)
        uNext(3:N-2) = uINext(3:end-3) * alpha + uINext(4:end-2) * (1-alpha);
        u(3:N-2) = uI(3:end-3) * alpha + uI(4:end-2) * (1-alpha);
    end
    
    if drawString == false && mod(t,1) == 0 && t >= start - 100
        lowRange = 28;
        highRange = 33;
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

winSize = 2^12;
overlap = 2^8;
% spectrogram(out*100,hanning(winSize),overlap,[],fs, 'MinThresh', -100,'yaxis')
plot(out);

% For loop
% if vectorRepresentation == false
%     for l = 2 : N - 1
%         if l > 2 && l < N - 1
%             % Newton Raphson
%             if l == floor(bowPoint)
% %                     if l == bowPoint
%                 alpha = (bowPoint - floor(bowPoint));
% %                     J = -1/h * (...
% %                         (alpha * (alpha - 1) * (alpha - 2))/ -6 * u(l-1) ...
% %                         + ((alpha - 1) * (alpha + 1) * (alpha - 2)) / 2 * u(l) ...
% %                         + (alpha * (alpha + 1) * (alpha - 2))/ -2 * u(l+1) ...
% %                         + (alpha * (alpha + 1) * (alpha - 1)) / 6 * u(l+2));
% 
%                 uI(2:end-1) = u(1:end-1) * (1-alpha) + u(2:end) * alpha;
%                 uIPrev(2:end-1) = uPrev(1:end-1) * (1-alpha) + uPrev(2:end) * alpha;
%                 bp = ceil(bowPoint);
% 
%                 b = 2/k * Vb - (2/k^2)*(uI(bp)-uIPrev(bp)) -(gamma/h)^2*(uI(bp+1) - 2*uI(bp) + uI(bp-1))...
%                     + kappa^2/h^4 * (uI(bp+2) - 4*uI(bp+1) + 6*uI(bp) - 4*uI(bp-1) + uI(bp-2)) + 2 * s0 * Vb ...
%                     - (2*s1/(k*h^2)) * ((uI(bp+1) - 2 * uI(bp) + uI(bp-1)) ...
%                       - (uIPrev(bp+1) - 2 * uIPrev(bp) + uIPrev(bp-1)));
%                 eps = 1;
%                 i = 0;
%                 while eps>tol
%                     q=qPrev-(1/h * Fb*A*qPrev*exp(-a*qPrev^2)+2*qPrev/k+b)/...
%                      (1/h * Fb*A*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k);
%                     eps = abs(q-qPrev);
%                     qPrev = q;
%                     i = i + 1;
%                     if i > 10000
%                         disp('whut')
%                     end
%                 end
%                 qSave(t) = q;
%                 uNextI(bp) = 2 * k * (q + Vb) + uIPrev(bp);
%             elseif l ~= ceil(bowPoint)
%                 uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
%                       - muSq * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) ...
%                       + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
%                       * ((u(l+1) - 2 * u(l) + u(l-1)) ...
%                       - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
%                       / (1 + s0 * k);
%             end
%         end
%     end