clear all;
close all;
clc;

fs = 44100;     % Sampling rate

E = 2E11;       % Young's Modulus
L = 0.40;       % String length
rho = 7850;     % Density of steel [kg/m^3]
T = 70;         % Tension [N]
c = 100; %sqrt(T / rho);% Wave speed
r = 0.001;      % String radius
A = r^2 * pi; % Cross-sectional area
m = A * L * rho; % Mass
inertia = pi / 2 * r^4;         % Moment of inertia
kappa = sqrt((E * inertia) / (rho * A * L^4));  % Stiffness

k = 1/fs;
h = c * k;

lambdaSq = (c * k / h)^2;

P = 1/3;

n = floor(L / h) + 1;

PIdx = floor (P*n);
uPrev = [0:1/PIdx:1 1-1/(n-PIdx):-1/(n-PIdx):0]';
u = uPrev;
plot (u);
n = n + 1;

uNext = zeros(n, 1);
Dxx = (sparse(2:n, 1:n-1, ones(1, n-1), n, n) + ...
    sparse(1:n, 1:n, -2 * ones(1, n), n, n) + ...
    sparse(1:n-1, 2:n, ones(1, n-1), n, n));
I= sparse(1:n, 1:n, ones(1, n), n, n);
B = 2 * I + lambdaSq * Dxx;

lengthSound = fs;
out = zeros(lengthSound, 1);
for t = 1 : lengthSound
    uNext = B * u - uPrev;
    out(t) = uNext(PIdx);
    uPrev = u;
    u = uNext;
    plot(uNext);
    ylim([-1 1]);
    drawnow;
end
plot(out);
% for t = 1 : fs / 10
%     for x = 2 : length(u) - 1
%         uNext(x) = 2 * u(x) - uPrev(x) + lambdaSq * ...
%             (u(x+1) - 2 * u(x) + u(x-1));
%     end
%     uPrev = u;
%     u = uNext;
%     if mod (t, 10) == 0
%         plot(uNext);
%         ylim([-1 1]);
%         drawnow;
%     end
% end