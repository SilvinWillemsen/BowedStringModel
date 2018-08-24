clear all;
close all;
clc;

fs = 44100;     % Sampling rate
E = 2E11;       % Young's Modulus
L = 0.64;       % String length
rho = 7850;     % Density of steel [kg/m^3]

r = 0.001;      % String radius
inertia = pi / 4 * r^4;         % Moment of inertia
kappa = sqrt((E * inertia) / rho);  % Stiffness

k = 1/fs;
h = sqrt(2 * kappa * k) * 1.2;
muSq = (k * kappa / h^2)^2;

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
% Dxxxx = (sparse(3:n, 1:n-2, ones(1, n-2), n, n) + ...
%         sparse(2:n, 1:n-1, -4 * ones(1, n-1), n, n) + ...
%         sparse(1:n, 1:n, 6 * ones(1, n), n, n) + ...
%         sparse(1:n-1, 2:n, -4 * ones(1, n-1), n, n) + ...
%         sparse(1:n-2, 3:n, ones(1, n-2), n, n));
Dxx = (sparse(2:n, 1:n-1, ones(1, n-1), n, n) + ...
    sparse(1:n, 1:n, -2 * ones(1, n), n, n) + ...
    sparse(1:n-1, 2:n, ones(1, n-1), n, n));
% Dxx(:,1:2) = 0;
% Dxx(:,end-1:end) = 0;
% Dxx(1:2,:) = 0;
% Dxx(end-1:end,:) = 0;
Dxxxx = Dxx * Dxx;

I = sparse(1:n, 1:n, ones(1, n), n, n);
B = 2 * I - muSq * Dxxxx;
Buse = B(3:n-2, 3:n-2);
lengthSound = fs * 10;
out = zeros(lengthSound, 1);
for t = 1 : lengthSound
    uNext(3:n-2) = B(3:n-2,3:n-2) * u(3:n-2) - uPrev(3:n-2);
    out(t) = uNext(length(uNext) - PIdx);
    uPrev = u;
    u = uNext;
    if mod(t,100) == 1
        plot(uNext);
        ylim([-1 1]);
        drawnow;
    end
end
plot(out);

for t = 1 : fs * 3
    for x = 3 : length(u) - 2
        uNext(x) = (2 - 6 * muSq) * u(x) + ...
                    4 * muSq * (u(x+1) + u(x-1)) - ...
                    muSq * (u(x-2) + u(x+2)) - uPrev(x);
    end
    out(t) = uNext(PIdx);
    uPrev = u;
    u = uNext;
end
plot(out);