clear all;
close all;
clc;

fs = 8000;     % Sampling rate
E = 2E11;       % Young's Modulus
L = 0.40;       % String length
rho = 7850;     % Density of steel [kg/m^3]
c = 21;        % Speed
gamma = c / L;
k = 1/fs;
r = 0.0001;      % String radius
inertia = pi / 4 * r^4;         % Moment of inertia
A = r^2 * pi; % Cross-sectional area
kappa = sqrt((E * inertia) / rho * A * L^4);  % Stiffness

% h = nthroot(k^2 * kappa^2 + 16 * gamma^2 * k^2, 4);
h = sqrt((gamma^2*k^2+sqrt(gamma^4 * k^4 + 16 * kappa^2 * k^2)) / 2);
muSq = (k * kappa / h^2)^2;
lambdaSq = (gamma * k / h);
% lambdaSq + 4 * muSq = 1;

P = 1/3;

n = floor(L / h) + 1;

cosWidth = round(n / 30);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * n);
u = zeros(n, 1);
u(floor(n * P - cosWidth / 2 : ceil(n * P + cosWidth / 2))) = raisedCos;
uPrev = u;
plot (u);
uNext = zeros(n,1);
Dxx = (sparse(2:n, 1:n-1, ones(1, n-1), n, n) + ...
    sparse(1:n, 1:n, -2 * ones(1, n), n, n) + ...
    sparse(1:n-1, 2:n, ones(1, n-1), n, n));
Dxxxx = Dxx * Dxx;

I = sparse(1:n, 1:n, ones(1, n), n, n);
B = 2 * I + lambdaSq * Dxx - muSq * Dxxxx;

Buse = B(3:n-2, 3:n-2);
lengthSound = fs;
out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

range = 3:n-2;
drawString = false;

for t = 1 : lengthSound
    uNext(range) = B(range, range) * u(range) - uPrev(range);
    out(t) = uNext(length(uNext) - PIdx);
    kinEnergy(t) = 1 / 2 * sum (h * ((1 / k * (u(range) - uPrev(range))).^2));
    potEnergy(t) = c^2 / 2 * sum (1 / h * ...
        (u(range + 1) - u(range)) .* (uPrev(range + 1) - uPrev(range)) + ...
        kappa^2 / 2 *  1 / h^3 * ...
        (u(range + 1) - 2 * u(range) + u(range - 1)) .* ...
        (uPrev(range + 1) - 2 * uPrev(range) + uPrev(range - 1)));
    
    uPrev = u;
    u = uNext;
    if mod(t,2) == 1 && drawString
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

for t = 1 : fs * 3
    for x = 3 : length(u) - 2
        uNext(x) = (2 - 2 * lambdaSq - 6 * muSq) * u(x) + ...
                   (lambdaSq + 4 * muSq) * (u(x+1) + u(x-1)) - ...
                    muSq * (u(x-2) + u(x+2)) - uPrev(x);
    end
    out(t) = uNext(PIdx);
    uPrev = u;
    u = uNext;
end
plot(out);