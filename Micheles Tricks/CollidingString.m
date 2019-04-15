clear all;
close all;
fs = 44100;
k = 1/fs;

drawThings = true;
drawSpeed = 1;


lengthSound = 2 * fs;

drawStart = 0;
f0 = 196.00;    % G3
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;     % Wave-Speed
T = c^2 * rho * A;
k = 1 / fs;     % Time step
E = 2e11;
I = r^4 * pi / 4;
L = 1;
kappa = sqrt(E*I/(rho*A));
s0 = 0;
s1 = 0;

h = sqrt((c^2*k^2 + 4*s1*k + sqrt((c^2*k^2 + 4*s1*k)^2+16*kappa^2*k^2)) / 2);
N = floor(1/h);
h = 1/N;
u = zeros(N, 1);
uPrev = zeros(N, 1);
%%  Raised cosine
width = 10;
loc = 0.25;
startIdx = floor(floor(loc * N) - width / 2);
endIdx = floor(floor(loc * N) + width / 2);

u(startIdx : endIdx) = u(startIdx : endIdx) + (1 - cos(2 * pi * [0:width]' / width)) / 2;
uPrev = u;

uNext = zeros(N, 1);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
colEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

alpha = 1.3;

% b = -10 * ones(N,1);
b = -1e-1-1e-1*(0:1/N:1-1/N)'-1*(0:1/N:1-1/N)'.^2;
% b(floor(N / 2) - 3: floor(N / 2)+ 3 ) = 0;

K = 5 * 10^8;

phiPrev = 0;
etaPrev = b - uPrev;
psiPrev = zeros(N, 1);
eta = b - uPrev;
etaNext = b - u;
Adiv = zeros(N, 1);
v = zeros(N, 1);

vec = 3:N-2;
eVec = 2:N-1; % energy vector
for n = 2:lengthSound
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / 2 * 1/h * sum((u(3:N) - u(2:N-1)) .* (uPrev(3:N) - uPrev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) ...
        .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)));
%     potEnergy2(n) = T / 2 * sum (1 / hS * (u(eVec + 1) - u(eVec)) .* (uPrev(eVec + 1) - uPrev(eVec))) ...
%         + (E * I) / 2 * 1/h^3 * sum((DS * u(1:N)) .* (DS * uPrev(1:N)));
    colEnergy(n) = 1/2 * sum(h * psiPrev.^2);
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + colEnergy(n);
    etaPrev = eta;
    eta = etaNext; 
    
    %% Calculate g
    if alpha == 1
        g = 0;
        if eta > 0
            g = ones(N,1) * sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta).^((alpha - 1)/2);
    end

    Adiv(vec) = rho * A / k^2 + g(vec).^2/4;
    v(vec) = rho * A / k^2 * (2 * u(vec) - uPrev(vec)) + T / h^2 * (u(vec+1) - 2 * u(vec) + u(vec-1)) ...
         - E * I / h^4 * (u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2))...
         + g(vec).^2/4 .* uPrev(vec) + psiPrev(vec) .* g(vec);
    uNext(vec) = v(vec) ./ Adiv(vec);
%     uNext(vec) = 2 * u(vec) - uPrev(vec) + T * k^2 / (rho * A * h^2) * (u(vec+1) - 2 * u(vec) + u(vec-1))...
%         - E * I * k^2 / (rho * A * h^4) * (u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2));
%       uNext = BS * u + CS * uPrev; 
    etaNext = b - uNext;
%     %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
    psiPrev = psi; 
    
    %% Update states
    uPrev = u;
    u = uNext;
    
    %% Draw functions
    if mod(n,drawSpeed) == 0 && drawThings == true
        subplot(3,1,1);
        cla
        plot(uNext);
        hold on;
        plot(b)
        subplot(3,1,2);
%         energyPlot = (totEnergy(10:n) / totEnergy(10)) - 1;
        cla
%         plot(kinEnergy(10:n) / kinEnergy(10) - 1);
%         hold on;
%         plot(potEnergy(10:n) / potEnergy(10) - 1);
        plot(totEnergy(10:n) / totEnergy(10) - 1);
        subplot(3,1,3);
        plot(eta)
       
        drawnow;
    end
end
plot(uNext)