clear all;
close all;

fs = 44100;
k = 1/fs;

drawThings = true;
drawSpeed = 5;
collision = true;
lengthSound = 1 * fs;
drawStart = lengthSound - 2 * drawSpeed;
% drawStart = 0;

f0 = 110.00;    
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;     % Wave-Speed
T = c^2 * rho * A;
k = 1 / fs;     % Time step
E = 2e11;
I = r^4 * pi / 4;
L = 1;
kappa = sqrt (E*I / (rho*A));
s0 = 1;
s1 = 0.00;
scaleFac = rho * A;
% scaleFac = 1;
h = sqrt((c^2*k^2 + 4*s1*scaleFac*k + sqrt((c^2*k^2 + 4*s1*scaleFac*k)^2+16*kappa^2*k^2)) / 2);
N = floor(1/h);
h = 1/N;

u = zeros(N, 1);
uPrev = zeros(N, 1);
%%  Raised cosine
width = 20;
loc = 1/pi;
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
rOCdamp0Energy = zeros(lengthSound, 1);
rOCdamp1Energy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);

alpha = 1.3;

b = -1 * ones(N,1);
% b = -1e-1-1e-1*(0:1/N:1-1/N)'-1*(0:1/N:1-1/N)'.^2;
width = floor(N / 6);
loc = floor(N / 5);
sharpness = 1;
closeness = 0.1;
if collision
    b(loc - floor(width / 2): loc + ceil(width / 2) ) = ((1 - cos(pi * [0:width] / width + 0.5 * pi)) * 0.5) * sharpness - sharpness - closeness;
end
K = 5 * 10^15;
% plot(b)
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
    %% Calculate energy of the system
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / 2 * 1/h * sum((u(3:N) - u(2:N-1)) .* (uPrev(3:N) - uPrev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) ...
        .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)));
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

    %% Update FDS
    Adiv(vec) = rho * A / k^2 + g(vec).^2/4 + scaleFac * s0/k;
    v(vec) = rho * A / k^2 * (2 * u(vec) - uPrev(vec)) + T / h^2 * (u(vec+1) - 2 * u(vec) + u(vec-1)) ...
         - E * I / h^4 * (u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2))...
         + scaleFac * s0 / k * uPrev(vec) ...
         + scaleFac * 2 * s1 / (h^2 * k) * (u(vec+1) - 2 * u(vec) + u(vec-1)) .* (uPrev(vec+1) - 2 * uPrev(vec) + uPrev(vec-1))...
         + g(vec).^2/4 .* uPrev(vec) + psiPrev(vec) .* g(vec);
    uNext(vec) = v(vec) ./ Adiv(vec);
    etaNext = b - uNext;
    
    %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
    
    %% Calculate rate-of-changes in energies (inner product of delta tdot with scheme) 
    rOCkinEnergy(n) = h * rho * A / (2 * k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = h * T / (2*k*h^2) * sum((u(vec+1) - 2 * u(vec) + u(vec-1)).* (uNext(vec) - uPrev(vec))) ...
         - h * E * I / (2 * k * h^4) * sum((u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2)) .* (uNext(vec) - uPrev(vec)));%...
    rOCdamp0Energy(n) = -scaleFac * 2 * s0 * h / (4 * k^2) * sum((uNext - uPrev).^2);
    rOCdamp1Energy(n) = scaleFac * 2 * h * s1 / (2 * k^2 * h^2) * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)) .* (uNext(eVec) - uPrev(eVec)));
    rOCcolEnergy(n) = h / (4 * k) * sum(g .* (psi + psiPrev) .* (uNext - uPrev));
    rOCTotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCdamp0Energy(n) - rOCdamp1Energy(n) - rOCcolEnergy(n); %including damping so should be 0
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
        plot(rOCdamp0Energy(10:n));
        plot(rOCdamp1Energy(10:n));
        plot(rOCcolEnergy(10:n));
        title("Rate-of-Change")
        legend(["Kin", "Pot", "$s_0$", "$s_1$", "Coll"], 'interpreter', 'latex')
        subplot(2,2,4)
        plot(rOCTotEnergy(10:n));
%         hold on;
%         plot(rOCdampEnergy(10:n));
        title("Total Rate-of-Change in energy");
        drawnow
    end
end
if ~drawThings
    plot(out)
end