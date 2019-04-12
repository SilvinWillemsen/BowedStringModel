clear all;
close all;
fs = 44100;
k = 1/fs;

f0 = 110;
w0 = 2 * pi * f0;
M = 0.01;
% M = k * K / ((2 * pi * f0)^2);

lengthSound = 2 * fs;

u = -1;
uPrev = -1;

%% Simple mass-spring
% for n = 2:lengthSound
%     uNext(n) = 2 * u - uPrev - k^2 * (K/M)^2 * u;
%     uPrev = u;
%     u = uNext(n);
% end
alpha = 1.3;
K = 5 * 10^6;
b = 0;
phiPrev = 0;
etaPrev = uPrev - b;
psiPrev = 0;
eta = 0;
for n = 2:lengthSound
    energy(n) = M * (1/k * (u - uPrev))^2 / 2 + M * w0^2*u*uPrev / 2 + psiPrev^2 / 2;
    etaPrev = eta;
    eta = u - b;
%     b = 0.1 * sin(2 * pi * n / fs);
    bSave(n) = b;
    %% Calculate psi dot / eta dot
    if alpha == 1
        g = 0;
        if eta > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end

    gSave(n) = g;
    %% Update the scheme
    A = M / k^2 + g^2 / 4;
    v = M / k^2 * (2 * u - uPrev) - M * w0^2 * u + g^2 / 4 * uPrev - psiPrev * g;
    uNext(n) = v / A;
  
    
    %% Update Psi
    psi = psiPrev + 0.5 * g * ((uNext(n) - b) - etaPrev);
    psiPrev = psi; 
    
    %% Update states
    uPrev = u;
    u = uNext(n);
    
    %% Draw functions
    if mod(n,100) == 0
        clf
        subplot(2,1,1);
        plot(uNext(1:n));
        hold on; plot(bSave(1:n));
        subplot(2,1,2);
        plot(energy(1:n));
        drawnow
    end
end
plot(uNext)