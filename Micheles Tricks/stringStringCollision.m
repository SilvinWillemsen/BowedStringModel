clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = true;
drawSpeed = 1;
lengthSound = 1 * fs;
drawStart = 0;

%% String 1 Variables
f1 = 110.00;    
rho1 = 7850;
r1 = 0.001;
A1 = r1^2 * pi;
c1 = f1 * 2;         % Wave-Speed
T1 = c1^2 * rho1 * A1;  % Tension
E1 = 0;           % Young's modulus
I1 = r1^4 * pi / 4;   % Area moment of inertia
L1 = 1;              % String Length
kappa1 = sqrt (E1*I1 / (rho1*A1));   % Stiffness coefficient

% Damping coefficients
s01 = rho1 * A1 * 0;
s11 = 0;
scaleFac1 = 1; % (scaling with mass/unit length)

% Grid spacing
h1 = sqrt((c1^2*k^2 + 4*s11*scaleFac1*k + sqrt((c1^2*k^2 + 4*s11*scaleFac1*k)^2+16*kappa1^2*k^2)) / 2);
N1 = floor(1/h1);
h1 = 1/N1;

u1 = zeros(N1, 1);
u1Prev = zeros(N1, 1);
u1Next = zeros(N1, 1);

%% String 1 Variables
f2 = 110.00;    
rho2 = 7850;
r2 = 0.001;
A2 = r2^2 * pi;
c2 = f2 * 2;         % Wave-Speed
T2 = c2^2 * rho2 * A2;  % Tension
E2 = 0;           % Young's modulus
I2 = r2^4 * pi / 4;   % Area moment of inertia
L2 = 1;              % String Length
kappa2 = sqrt (E2*I2 / (rho2*A2));   % Stiffness coefficient

% Damping coefficients
s02 = rho2 * A2 * 0;
s12 = 0;
scaleFac2 = 1; % (scaling with mass/unit length)

% Grid spacing
h2 = sqrt((c2^2*k^2 + 4*s12*scaleFac2*k + sqrt((c2^2*k^2 + 4*s12*scaleFac2*k)^2+16*kappa2^2*k^2)) / 2);
N2 = floor(1/h2);
h2 = 1/N2;

offset = 0.1;

u2 = zeros(N2, 1) - offset;
u2Prev = zeros(N2, 1) - offset;
u2Next = zeros(N2, 1) - offset;

%% Collision Variables
alpha = 1.3;
K = 5 * 10^6;
cL = floor (N1 / 2); % collision location

%% Initialise
etaPrev = u2 - u1;
psiPrev = zeros(N1, 1);
psi = zeros(N1, 1);
eta = u2Prev - u1Prev;
etaNext = u2 - u1;

%% Excitation (raised cosine)
width = 20;
loc = 1/2;
startIdx = floor(floor(loc * N1) - width / 2);
endIdx = floor(floor(loc * N1) + width / 2);
u1(startIdx : endIdx) = u1(startIdx : endIdx) + (1 - cos(2 * pi * [0:width]' / width)) / 2;
u1Prev = u1;

%% Initialise Energy Vectors
kinEnergy1 = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);
kinEnergy2 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);
colEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

rOCkinEnergy1 = zeros(lengthSound, 1);
rOCpotEnergy1 = zeros(lengthSound, 1);
rOCcolEnergy1 = zeros(lengthSound, 1);
rOCdamp0Energy1 = zeros(lengthSound, 1);
rOCdamp1Energy1 = zeros(lengthSound, 1);
rOCtotEnergy1 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCcolEnergy2 = zeros(lengthSound, 1);
rOCdamp0Energy2 = zeros(lengthSound, 1);
rOCdamp1Energy2 = zeros(lengthSound, 1);
rOCtotEnergy2 = zeros(lengthSound, 1);

rOCtotEnergy = zeros(lengthSound, 1);

A = zeros(N1, 1);
v = zeros(N1, 1);

vec1 = 3:N1-2;
eVec1 = 2:N1-1; % energy vector
vec2 = 3:N2-2;
eVec2 = 2:N2-1; % energy vector

etaNextDiv = zeros(N1 ,1);
etaNextV = zeros(N1 ,1);

outputPos = floor(N1 / 3);
scaleH = 1;
for n = 2:lengthSound
    %% Calculate energy of the system
    kinEnergy1(n) = rho1 * A1 / 2 * h1 * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T1 / 2 * 1/h1 * sum((u1(3:N1) - u1(2:N1-1)) .* (u1Prev(3:N1) - u1Prev(2:N1-1)))...
        + E1 * I1 / 2 * 1/h1^3 * sum((u1(eVec1+1) - 2 * u1(eVec1) + u1(eVec1-1)) ...
        .* (u1Prev(eVec1+1) - 2 * u1Prev(eVec1) + u1Prev(eVec1-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    kinEnergy2(n) = rho2 * A2 / 2 * h2 * sum((1/k * (u2 - u2Prev)).^2);
    potEnergy2(n) = T2 / 2 * 1/h2 * sum((u2(3:N2) - u2(2:N2-1)) .* (u2Prev(3:N2) - u2Prev(2:N2-1)))...
        + E2 * I2 / 2 * 1/h2^3 * sum((u2(eVec2+1) - 2 * u2(eVec2) + u2(eVec2-1)) ...
        .* (u2Prev(eVec2+1) - 2 * u2Prev(eVec2) + u2Prev(eVec2-1)));
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    colEnergy(n) = 1/2 * h1 * sum(psiPrev.^2);
    totEnergy(n) = energy1(n) + energy2(n) + colEnergy(n);
    
    etaPrev = eta;
    eta = etaNext; 
    
    %% Calculate g
    if alpha == 1
        g = 0;
        if eta > 0
            g = ones(N1,1) * sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta).^((alpha - 1)/2);
    end

    %% Calculate etaNext (STILL NEED TO INCLUDE DAMPING)
    etaNextDiv(vec1) =  1 + g(vec2).^2 / 4 * 1 / ((rho2 * A2) / k^2  + s02 / k) + g(vec1).^2 / 4 * 1 / ((rho1 * A1) / k^2 + s01 / k);
    etaNextV(vec1) = 2 * u2(vec2) - u2Prev(vec2) + (T2 / h2^2 * (u2(vec2+1) - 2 * u2(vec2) + u2(vec2-1))...
        - E2 * I2 / h2^4 * (u2(vec2+2) - 4 * u2(vec2+1) + 6 * u2(vec2) - 4 * u2(vec2-1) + u2(vec2-2))...
        + s02 / k * u2Prev(vec2)...
        + (g(vec2).^2 / 4 .* etaPrev(vec2) - psiPrev(vec2) .* g(vec2))) * k^2 / (rho2 * A2)...
        ...
        - (2 * u1(vec1) - u1Prev(vec1) + (T1 / h1^2 * (u1(vec1+1) - 2 * u1(vec1) + u1(vec1-1))...
        - E1 * I1 / h1^4 * (u1(vec1+2) - 4 * u1(vec1+1) + 6 * u1(vec1) - 4 * u1(vec1-1) + u1(vec1-2))...
        + s01 / k * u1Prev(vec1)...
        - (g(vec1).^2 / 4 .* etaPrev(vec1) - psiPrev(vec1) .* g(vec1))) * k^2 / (rho1 * A1));
    etaNext(vec1) = etaNextV(vec1) ./ etaNextDiv(vec1);
    
    %% Update FDS
%     Adiv(vec) = rho * A / k^2;
    u1Next(vec1) = (rho1 * A1 / k^2 * (2 * u1(vec1) - u1Prev(vec1)) + T1 / h1^2 * (u1(vec1+1) - 2 * u1(vec1) + u1(vec1-1)) ...
         - E1 * I1 / h1^4 * (u1(vec1+2) - 4 * u1(vec1+1) + 6 * u1(vec1) - 4 * u1(vec1-1) + u1(vec1-2))...
         + s01 / k  * u1Prev(vec1) ...
         + (g(vec1).^2/4 .* (etaNext(vec1) - etaPrev(vec1)) + psiPrev(vec1) .* g(vec1))) * 1 / ((rho1 * A1) / k^2 + s01 / k);
%               
    u2Next(vec2) = (rho2 * A2 / k^2 * (2 * u2(vec2) - u2Prev(vec2)) + T2 / h2^2 * (u2(vec2+1) - 2 * u2(vec2) + u2(vec2-1)) ...
         - E2 * I2 / h2^4 * (u2(vec2+2) - 4 * u2(vec2+1) + 6 * u2(vec2) - 4 * u2(vec2-1) + u2(vec2-2))...
         + s02 / k  * u2Prev(vec2) ...
         - (g(vec2).^2/4 .* (etaNext(vec2) - etaPrev(vec2)) + psiPrev(vec2) .* g(vec2))) * 1 / ((rho2 * A2) / k^2 + s02 / k);

    %% Update Psi
    psi(vec1) = psiPrev(vec1) + 0.5 * g(vec1) .* (etaNext(vec1) - etaPrev(vec1));
    
%     %% Calculate rate-of-changes in energies for checking damping (inner product of delta tdot with scheme) 
    rOCkinEnergy1(n) = h1 * rho1 * A1 / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = h1 * T1 / (2*k*h1^2) * sum((u1(vec1+1) - 2 * u1(vec1) + u1(vec1-1)).* (u1Next(vec1) - u1Prev(vec1))) ...
         - h1 * E1 * I1 / (2 * k * h1^4) * sum((u1(vec1+2) - 4 * u1(vec1+1) + 6 * u1(vec1) - 4 * u1(vec1-1) + u1(vec1-2)) .* (u1Next(vec1) - u1Prev(vec1)));%...
    rOCdamp0Energy1(n) = -scaleFac1 * 2 * s01 * h1 / (4 * k^2) * sum((u1Next - u1Prev).^2);
    rOCdamp1Energy1(n) = scaleFac1 * 2 * h1 * s11 / (2 * k^2 * h1^2) * sum((u1(eVec1+1) - 2 * u1(eVec1) + u1(eVec1-1)) .* (u1Prev(eVec1+1) - 2 * u1Prev(eVec1) + u1Prev(eVec1-1)) .* (u1Next(eVec1) - u1Prev(eVec1)));
    rOCcolEnergy1(n) = h1 / (4 * k) * sum(g .* (psi + psiPrev) .* (u1Next - u1Prev));
    rOCTotEnergy1(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0Energy1(n) - rOCdamp1Energy1(n) - rOCcolEnergy1(n);
    
    rOCkinEnergy2(n) = h2 * rho2 * A2 / (2 * k^3) * sum((u2Next - 2 * u2 + u2Prev) .* (u2Next - u2Prev));
    rOCpotEnergy2(n) = h2 * T2 / (2*k*h2^2) * sum((u2(vec2+1) - 2 * u2(vec2) + u2(vec2-1)).* (u2Next(vec2) - u2Prev(vec2))) ...
         - h2 * E2 * I2 / (2 * k * h2^4) * sum((u2(vec2+2) - 4 * u2(vec2+1) + 6 * u2(vec2) - 4 * u2(vec2-1) + u2(vec2-2)) .* (u2Next(vec2) - u2Prev(vec2)));%...
    rOCdamp0Energy2(n) = -scaleFac2 * 2 * s02 * h2 / (4 * k^2) * sum((u2Next - u2Prev).^2);
    rOCdamp1Energy2(n) = scaleFac2 * 2 * h2 * s12 / (2 * k^2 * h2^2) * sum((u2(eVec2+1) - 2 * u2(eVec2) + u2(eVec2-1)) .* (u2Prev(eVec2+1) - 2 * u2Prev(eVec2) + u2Prev(eVec2-1)) .* (u2Next(eVec2) - u2Prev(eVec2)));
    rOCcolEnergy2(n) = -h2 / (4 * k) * sum(g .* (psi + psiPrev) .* (u2Next - u2Prev));
    rOCTotEnergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n) - rOCdamp0Energy2(n) - rOCdamp1Energy2(n) - rOCcolEnergy2(n);
    
    rOCTotEnergy(n) = rOCTotEnergy1(n) + rOCTotEnergy2(n); %including damping so should be 0
    %% Update states
    psiPrev = psi; 
    u1Prev = u1;
    u1 = u1Next;
    u2Prev = u2;
    u2 = u2Next;
    
    out1(n) = u1Next(outputPos);
    out2(n) = u2Next(outputPos);
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(2,1,1);
        cla
        plot(u1Next);
        hold on;
        plot(u2Next);
       
        subplot(2,1,2);
        cla
%         plot(totEnergy(10:n) / totEnergy(10) - 1)
plot(energy1(10:n));
hold on;
plot(energy2(10:n));
%         plot(rOCTotEnergy(10:n));
% %         plot(psi)
%         hold on; 
%         plot(rOCdamp0Energy1(10:n))
%         plot(rOCdamp0Energy2(10:n))
%         plot((rOCcolEnergy1(10:n)+rOCcolEnergy2(10:n)) / 2);
%         plot(energy1(10:n) - energy1(10));
%         plot(energy2(10:n) - energy2(10));
%         plot(colEnergy(10:n) - colEnergy(10));
%         plot(colEnergy(100:n) / colEnergy(100) - 1)
%         if s0 == 0 && s1 == 0
%             plot(totEnergy(10:n) / totEnergy(10) - 1);
%             title("Total normalised energy");
%         else
%             plot(totEnergy(10:n));
%             title("Total energy")
%         end
%         subplot(2,2,3);
%         cla
% %         totEnergyPlot(n) = totEnergy(n) / totEnergy(10) - 1;
%         plot(rOCkinEnergy(10:n));
%         hold on;
%         plot(rOCpotEnergy(10:n), '--');
%         plot(rOCdamp0Energy(10:n));
%         plot(rOCdamp1Energy(10:n));
%         plot(rOCcolEnergy(10:n));
%         title("Rate-of-Change")
%         legend(["Kin", "Pot", "$s_0$", "$s_1$", "Coll"], 'interpreter', 'latex')
%         subplot(2,2,4)
%         plot(rOCTotEnergy(10:n));
% %         hold on;
% %         plot(rOCdampEnergy(10:n));
%         title("Total Rate-of-Change in energy");
        drawnow
    end
end
if ~drawThings
    plot(out1)
    hold on;
    plot(out2 + (max(out1) - min(out1)))
end