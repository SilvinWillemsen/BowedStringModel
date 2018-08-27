clear all;
close all;
clc;

fs = 8000;      % Sampling rate
k = 1 / fs;     % Time step
L = 1;          % String length
rho = 7850;     % Density of steel [kg/m^3]
r = 0.001;      % String radius
A = pi * r^2;   % Cross-sectional area
E = 2E11;       % Young's Modulus
inertia = pi / 4 * r^4; % Moment of inertia
kappa = sqrt((E * inertia) / (rho * A * L^4));  % Stiffness

h = sqrt(2 * kappa * k); % Grid spacing
muSq = (k * kappa / h^2)^2; % Courant number squared

P = 1/2; % plucking position

n = floor(L / h) + 1; % Number of grid-points

%% Raised cosine input
cosWidth = round(n / 20);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * n);

%% Initialise state vectors
u = zeros(n, 1);
u(floor(n * P - cosWidth / 2 : ceil(n * P + cosWidth / 2))) = raisedCos;
uPrev = u;
uNext = zeros(n,1);

%% Matrix representation

Dxx = (sparse(2:n, 1:n-1, ones(1, n-1), n, n) + ...
    sparse(1:n, 1:n, -2 * ones(1, n), n, n) + ...
    sparse(1:n-1, 2:n, ones(1, n-1), n, n));
Dxxxx = Dxx * Dxx;

I = sparse(1:n, 1:n, ones(1, n), n, n);
B = 2 * I - muSq * Dxxxx;
Buse = B(3:n-2, 3:n-2);
lengthSound = fs / 10;
out = zeros(lengthSound, 1);

drawBar = false;
range = 3 : n-2;

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

for t = 1 : lengthSound
    uNext(range) = Buse * u(range) - uPrev(range);
    out(t) = uNext(length(uNext) - PIdx);
    kinEnergy(t) = 1 / 2 * sum (h * ((1 / k * (u(range) - uPrev(range))).^2));
    potEnergy(t) = kappa^2 / 2 * sum (1 / h^3 * ...
        (u(range + 1) - 2 * u(range) + u(range - 1)) .* ...
        (uPrev(range + 1) - 2 * uPrev(range) + uPrev(range - 1)));
    uPrev = u;
    u = uNext;
    
    if mod(t,2) == 1 && drawBar
        plot(uNext);
        ylim([-1 1]);
        drawnow;
    end
end
plot(kinEnergy);
hold on; plot(potEnergy);
totEnergy = kinEnergy + potEnergy;
plot(totEnergy);
plot(out);

%% For-loop representation
% for t = 1 : fs * 3
%     for x = 3 : length(u) - 2
%         uNext(x) = (2 - 6 * muSq) * u(x) + ...
%                     4 * muSq * (u(x+1) + u(x-1)) - ...
%                     muSq * (u(x-2) + u(x+2)) - uPrev(x);
%     end
%     out(t) = uNext(PIdx);
%     uPrev = u;
%     u = uNext;
% end
% plot(out);