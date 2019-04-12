clear all;
close all;
fs = 44100;
k = 1/fs;

drawThings = false;
drawSpeed = 100;

f1 = 4;
f2 = 100;
w1 = 2 * pi * f1;
w2 = 2 * pi * f2;
M1 = 1;
M2 = 0.0001;

lengthSound = 2 * fs;

u1 = 1;
u1Prev = 1;
u2 = -1;
u2Prev = -1;

%% Initialise vectors for speed
u1Next = zeros(lengthSound, 1);
u2Next = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

alpha = 1.3;
K = 5 * 10^6;

phiPrev = 0;
etaPrev = u2 - u1;
psiPrev = 0;
eta = u2Prev - u1Prev;
etaNext = u2 - u1;

for n = 2:lengthSound
    energy1(n) = M1 * (1/k * (u1 - u1Prev))^2 / 2 + M1 * w1^2*u1*u1Prev / 2;
    energy2(n) = M2 * (1/k * (u2 - u2Prev))^2 / 2 + M2 * w2^2*u2*u2Prev / 2;
    totEnergy(n) = energy1(n) + energy2(n) + psiPrev^2 / 2;
    
    etaPrev = eta;
    eta = etaNext; 

    %% Calculate g
    if alpha == 1
        g = 0;
        if eta > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end
    
    %% Update etaNext
    etaNext = (2 * u2 - u2Prev + (-M2 * w2^2 * u2 + g^2/4 * etaPrev - psiPrev * g) * k^2 / M2 ...
        - (2 * u1 - u1Prev + (-M1 * w1^2 * u1 - g^2 / 4 * etaPrev + psiPrev * g) * k^2 / M1)) ...
        / (1 + g^2/4 * k^2 / M2 + g^2/4 * k^2/M1);

    u1Next(n) = 2 * u1 - u1Prev + (-M1 * w1^2 * u1 + (g^2/4 * (etaNext - etaPrev) + psiPrev * g)) * k^2/M1;
    u2Next(n) = 2 * u2 - u2Prev + (-M2 * w2^2 * u2 - (g^2/4 * (etaNext - etaPrev) + psiPrev * g)) * k^2/M2;
    
    %% Update Psi
    psi = psiPrev + 0.5 * g * (etaNext - etaPrev);
    psiPrev = psi; 
    
    %% Update states
    u1Prev = u1;
    u1 = u1Next(n);
    
    u2Prev = u2;
    u2 = u2Next(n);
    
    %% Draw functions
    if mod(n,drawSpeed) == 0 && drawThings == true
        clf
        subplot(2,1,1);
        plot(u1Next(1:n));
        hold on; plot(u2Next(1:n));
        subplot(2,1,2);
        totEnergyPlot = (totEnergy(10:n) / totEnergy(10)) - 1;
        plot(totEnergyPlot);
        drawnow
    end
end
plot(u2Next)