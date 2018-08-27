clear all;
close all;
clc;

fs = 8000;     % Sampling rate
k = 1 / fs;     % Time step
L = 0.40;       % String length
rho = 7850;     % Density of steel [kg/m^3]
c = 100;         % sqrt(T / rho);% Wave speed
h = c * k;      % Grid spacing

lambdaSq = (c * k / h)^2; % Courant number squared

P = 1/2;              % Plucking position

n = floor(L / h) + 1; % Number of grid-points

%% Raised cosine input
cosWidth = round(n / 10);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * n);
u = zeros(n, 1);
u(floor(n * P - cosWidth / 2 : ceil(n * P + cosWidth / 2))) = raisedCos;
uPrev = u;

%% Matrix Representation
uNext = zeros(n, 1);
Dxx = (sparse(2:n, 1:n-1, ones(1, n-1), n, n) + ...
    sparse(1:n, 1:n, -2 * ones(1, n), n, n) + ...
    sparse(1:n-1, 2:n, ones(1, n-1), n, n));
I = sparse(1:n, 1:n, ones(1, n), n, n);
B = 2 * I + lambdaSq * Dxx;

lengthSound = fs / 10;
out = zeros(lengthSound, 1);
kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);


range = 2 : n - 1;
for t = 1 : lengthSound
    uNext(range) = B(range,range) * u(range) - uPrev(range);
    kinEnergy(t) = sum (1 / 2 * h * ((1 / k * (u(range) - uPrev(range))).^2));
    potEnergy(t) = sum (c^2 / 2 * 1 / h * ...
        (u(range + 1) - u(range)) .* (uPrev(range + 1) - uPrev(range)));
    out(t) = uNext(PIdx);

    uPrev = u;
    u = uNext;
end
totEnergy = kinEnergy + potEnergy;
plot(totEnergy);
 
% %% For-loop representation
% potEnergy = zeros(lengthSound, 1);
% kinEnergy = zeros(lengthSound, 1);
% 
% for t = 1 : lengthSound
%     for l = 2 : length(u) - 1
%         uNext(l) = 2 * (1 - lambdaSq) * u(l) + lambdaSq * (u(l - 1) + u(l + 1)) - uPrev(l);
%         kinEnergy(t) = kinEnergy(t) + 1 / 2 * h * ((1/ k * (u(l) - uPrev(l)))^2);
%         potEnergy(t) = potEnergy(t) + c^2 / 2 * 1 / h * ...
%             (u(l+1) - u(l)) * (uPrev(l + 1) - uPrev(l));
%     end
%     uPrev = u;
%     u = uNext;
% 
%     out(t) = uNext(PIdx);
% end
% 
% totEnergy = kinEnergy + potEnergy;