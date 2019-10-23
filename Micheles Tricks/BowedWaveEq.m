clear all;
% close all;
figure
drawString = false;
drawSpeed = 5;

fs = 44100;
k = 1 / fs;     % Time step

%% Extra Settings
lengthSound = fs * 5;

%% String variables
L = 1;
f0 = 100;    
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rho * A;  % Tension
h = c * k;      % Grid spacing

N = floor(L / h); % Number of grid-points
h = L / N; % Recalculate gridspacing

lambdaSq = (c * k / h)^2 % Courant number squared

BP = 1/pi;      % Bowing Position

%% Initialise state vectors
u = zeros(N, 1);
uNext = zeros(N, 1);

%% Raised cosine
% amp = 0.1;
% width = 10;
% loc = 2/5;
% startIdx = floor(floor(loc * N) - width / 2);
% endIdx = floor(floor(loc * N) + width / 2);
% u(startIdx : endIdx) = u(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;

uPrev = u;

u2 = zeros(N, 1);
u2Next = zeros(N, 1);
%% Raised cosine
% amp = 0.1;
% width = 10;
% loc = 2/5;
% startIdx = floor(floor(loc * N) - width / 2);
% endIdx = floor(floor(loc * N) + width / 2);
% u2(startIdx : endIdx) = u(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;

u2Prev = u2;

%% Matrix Representation
Ndec = N;
N = N - 4;
Ndec = Ndec - N;
Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
    sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
    sparse(1:N-1, 2:N, ones(1, N-1), N, N));
I = sparse(1:N, 1:N, ones(1, N), N, N);

B = 2 * I + lambdaSq * Dxx;
kappa = 0;
kOh = (kappa/h^2)^2;
gOh = (c/h)^2;
s0 = 0;
s1 = 0;
phi = (2 * s1)/(k * h^2);
bB = (sparse(3:N, 1:N-2, kOh * ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (6*kOh + 2 * gOh - 2/k^2 + 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, kOh * ones(1, N-2), N, N));
bC = sparse(2:N, 1:N-1, phi * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (2/k^2 - 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, phi * ones(1, N-1), N, N);
N = N + 4;
matVec = 1 + Ndec / 2 : N - Ndec / 2;
BB = zeros(N);
BB(matVec, matVec) = B;
B = BB;

CC = zeros(N);
CC(matVec, matVec) = c;
C = CC;

bBB = zeros(N);
bBB(matVec, matVec) = bB;
bB = bBB;

bCC = zeros(N);
bCC(matVec, matVec) = bC;
bC = bCC;

VbInit = 0.2;
Fb = 0.01;
a = 100;

BM = sqrt(2 * a) * exp(1/2);
bp = floor(BP * N);

% qPrev = -Vb;
qPrev = -VbInit;
qPrev2 = -VbInit;
I = zeros(N,1);
I(bp) = 1;
rampVal = 10000;
out = zeros(lengthSound, 1);
J = 1/h * I;

tol = 1e-7;

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
energy = zeros(lengthSound, 1);

rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCbowEnergy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);

kinEnergy2 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCbowEnergy2 = zeros(lengthSound, 1);
rOCtotEnergy2 = zeros(lengthSound, 1);


for n = 1 : lengthSound
    if n < rampVal
        Vb = VbInit * n / rampVal;
    else
        Vb = VbInit;
    end
    
    % Newton-Raphson
    b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
    eps = 1;
    i = 0;
    while eps>tol
        q=qPrev-(1/h * Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0*qPrev+b)/...
         (1/h * Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0);
        eps = abs(q-qPrev);
        qPrev = q;
        i = i + 1;
        if i > 10000
            disp('Nope');
        end
    end
    qSave(n) = q;
    
%     q = 2 / k * (u(bp) - uPrev(bp)) - qPrev - 2 * Vb;
    %% simple backwards
%     q = (u2(bp) - uPrev2(bp)) / k - Vb;
    
    excitation = k^2 * J * Fb * sqrt(2 * a) * q * exp(-a*q^2 + 1/2);
    uNext = B * u - uPrev - excitation / (rho * A);

    eVec = 2:N-1;
    
    kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2*h) * sum((u(eVec + 1) - u(eVec)) .* (uPrev(eVec + 1) - uPrev(eVec)));
    energy(n) = kinEnergy(n) + potEnergy(n);
    
    rOCkinEnergy(n) = rho * A * h / (2 * k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = T / (2 * k * h) * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) .* (uNext(eVec) - uPrev(eVec)));
    rOCbowEnergy(n) = -Fb * (uNext(bp) - uPrev(bp)) / (2 * k) * sqrt(2*a) * q * exp(-a*q^2 + 1/2);
    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCbowEnergy(n);
    
    
    q2 = 2 / k * (u2(bp) - u2Prev(bp)) - qPrev2 - 2 * Vb;

    Bm = Fb * (sqrt(2 * a) * exp(-a * q2^2 + 1/2));
    u2Next = ((rho * A / k^2) * (B * u2 - u2Prev) + J * Bm .* (u2Prev / (2 * k) + Vb)) ./ (rho * A / k^2 + J * Bm / (2 * k));
%     u2Next = ((rho * A / k^2) * (B * u2 - u2Prev)) ./ (rho * A / k^2);
    
    kinEnergy2(n) = rho * A * h / 2 * sum((1/k * (u2 - u2Prev)).^2);
    potEnergy2(n) = T / (2*h) * sum((u2(eVec + 1) - u2(eVec)) .* (u2Prev(eVec + 1) - u2Prev(eVec)));
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    rOCkinEnergy2(n) = rho * A * h / (2 * k^3) * sum((u2Next - 2 * u2 + u2Prev) .* (u2Next - u2Prev));
    rOCpotEnergy2(n) = T / (2 * k * h) * sum((u2(eVec+1) - 2 * u2(eVec) + u2(eVec-1)) .* (u2Next(eVec) - u2Prev(eVec)));
    rOCbowEnergy2(n) = -Fb * sqrt(2*a) * ((u2Next(bp) - u2Prev(bp)) / (2 * k) - Vb) * exp(-a*q2^2 + 1/2) * (u2Next(bp) - u2Prev(bp)) / (2 * k) ;
    rOCtotEnergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n) - rOCbowEnergy2(n);
    
    if mod(n,drawSpeed) == 0 && drawString == true %&& n > 10000
        subplot(3,2,1);
        plot(uNext);
        
        subplot(3,2,2);
        plot(u2Next);
        
        subplot(3,2,3);
        plot(energy(10:n) / energy(10) - 1);
        
        subplot(3,2,4);
        plot(energy2(10:n) / energy2(10) - 1);
        
        subplot(3,2,5);
        plot(rOCtotEnergy(10:n))
        
        subplot(3,2,6);
        plot(rOCtotEnergy2(10:n))
        drawnow;
    end
    
    
    qPrev = q;
    uPrev = u;
    u = uNext;
    
    qPrev2 = q2;
    u2Prev = u2;
    u2 = u2Next;
    out(n) = u(floor(2*N / 3));
    out2(n) = u2(floor(2*N / 3));
end

% totEnergy = kinEnergy + potEnergy;
% totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
% maxTotEnergy = max(totEnergy) - min(totEnergy);
% % plot (totEnergy);
plot(out); 
hold on;
plot(out2)
legend(["Output NR","Output Tricks"])