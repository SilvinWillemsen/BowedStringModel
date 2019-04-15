clear all;
close all;

fs = 44100;
k = 1/fs;

drawThings = false;
drawSpeed = 100;
drawStart = 20000;

lengthSound = 2 * fs;

f0 = 110.00;    
rho = 7850;
r = 0.0001;
A = r^2 * pi;
c = f0 * 2;     % Wave-Speed
T = c^2 * rho * A;
k = 1 / fs;     % Time step
E = 0;
I = r^4 * pi / 4;
L = 1;
kappa = sqrt(E*I/(rho*A));
s0 = 1;
s1 = 0.00;

h = sqrt((c^2*k^2 + 4*s1*k + sqrt((c^2*k^2 + 4*s1*k)^2+16*kappa^2*k^2)) / 2);
N = floor(1/h);
h = 1/N;
% [BS, CS, NS, hS, DS, DS4] = newCreateString (c, kappa, L, s0, s1, k);

u = zeros(N, 1);
uPrev = zeros(N, 1);
%%  Raised cosine
width = 20;
loc = 0.45;
startIdx = floor(floor(loc * N) - width / 2);
endIdx = floor(floor(loc * N) + width / 2);

u(startIdx : endIdx) = u(startIdx : endIdx) + (1 - cos(2 * pi * [0:width]' / width)) / 2;
uPrev = u;

uNext = zeros(N, 1);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
colEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);
rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCcolEnergy = zeros(lengthSound, 1);
rOCdampEnergy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);

alpha = 1.3;

b = -10 * ones(N,1);
% b = -1e-1-1e-1*(0:1/N:1-1/N)'-1*(0:1/N:1-1/N)'.^2;

b(floor(N / 2) - floor(N / 4): floor(N / 2)+ floor(N / 4) ) = -0.1;

K = 5 * 10^6;

phiPrev = 0;
etaPrev = b - uPrev;
psiPrev = zeros(N, 1);
eta = b - uPrev;
etaNext = b - u;
Adiv = zeros(N, 1);
v = zeros(N, 1);

vec = 3:N-2;
eVec = 2:N-1; % energy vector

outputPos = floor(N / 3);
for n = 2:lengthSound
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / 2 * 1/h * sum((u(3:N) - u(2:N-1)) .* (uPrev(3:N) - uPrev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) ...
        .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)));
%     potEnergy2(n) = T / 2 * sum (1 / hS * (u(eVec + 1) - u(eVec)) .* (uPrev(eVec + 1) - uPrev(eVec))) ...
%         + (E * I) / 2 * 1/h^3 * sum((DS * u(1:N)) .* (DS * uPrev(1:N)));
%     energyLoss(n) = -2 * rho * A * (s0 * h * sum((1/(2*k) * (uNext - uPrev)).^2) + s1 * sum(1/k * (u(eVec + 1) - u (eVec)) .* (uPrev(eVec + 1) - uPrev(eVec))));
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

    Adiv(vec) = rho * A / k^2 + g(vec).^2/4 + rho * A * s0/k;
    v(vec) = rho * A / k^2 * (2 * u(vec) - uPrev(vec)) + T / h^2 * (u(vec+1) - 2 * u(vec) + u(vec-1)) ...
         - E * I / h^4 * (u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2))...
         + rho * A * s0 / k * uPrev(vec) ...
         + rho * A * 2 * s1 / (h^2 * k) * (u(vec+1) - 2 * u(vec) + u(vec-1)) .* (uPrev(vec+1) - 2 * uPrev(vec) + uPrev(vec-1))...
         + g(vec).^2/4 .* uPrev(vec) + psiPrev(vec) .* g(vec);
    uNext(vec) = v(vec) ./ Adiv(vec);
%     uNext(vec) = 2 * u(vec) - uPrev(vec) + T * k^2 / (rho * A * h^2) * (u(vec+1) - 2 * u(vec) + u(vec-1))...
%         - E * I * k^2 / (rho * A * h^4) * (u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2));
%       uNext = BS * u + CS * uPrev; 
    etaNext = b - uNext;
%     %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
    
    %% Calculate rate-of-changes in energies (inner product of delta tdot with scheme) 
    rOCkinEnergy(n) = h * rho * A / (2 * k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = h * T / (2*k*h^2) * sum((u(vec+1) - 2 * u(vec) + u(vec-1)).* (uNext(vec) - uPrev(vec))) ...
         - h * E * I / (2 * k * h^4) * sum((u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2)) .* (uNext(vec) - uPrev(vec)));%...
    rOCdampEnergy(n) = -rho * A * 2 * s0 * h / (4 * k^2) * sum((uNext - uPrev).^2);
    rOCcolEnergy(n) = h / (4 * k) * sum(g .* (psi + psiPrev) .* (uNext - uPrev));
    rOCTotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCdampEnergy(n) - rOCcolEnergy(n);
    %% Update states
    psiPrev = psi; 
    uPrev = u;
    u = uNext;
    out(n) = uNext(outputPos);
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(2,2,1);
        cla
        plot(uNext);
        hold on;
        plot(b)
        
        subplot(2,2,2);
        if s0 == 0 && s1 == 0
            plot(totEnergy(10:n) / totEnergy(10) - 1);
            title("Total normalised energy");
        else
            plot(totEnergy(10:n));
            title("Total energy")
        end
        subplot(2,2,3);
        cla
%         totEnergyPlot(n) = totEnergy(n) / totEnergy(10) - 1;
        plot(rOCkinEnergy(10:n));
        hold on;
        plot(rOCpotEnergy(10:n), '--');
        plot(rOCdampEnergy(10:n));
        plot(rOCcolEnergy(10:n));
        title("Rate-of-Change")
        legend(["Kin", "Pot", "Damp", "Coll"])
        subplot(2,2,4)
        plot(rOCTotEnergy(10:n));
        title("Total Rate-of-Change in energy");
        drawnow
    end
end
plot(out)