clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = true;
drawSpeed = 1000;
lengthSound = 5 * fs;
drawStart = 0;

%% String Variables
f0 = 200;    
rho = 7850;
r = 0.0009;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rho * A;  % Tension
k = 1 / fs;         % Time step
E = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappa = sqrt (E*I / (rho*A));   % Stiffness coefficient

% Damping coefficients
s0 = 0.0;
s1 = 0.0;

[B1, C1, N1, h1, Dxx, Dxxxx] = unscaledCreateString (rho, A, T, E, I, L, s0, s1, k);
[B2, C2, N2, h2, Dxx, Dxxxx] = unscaledCreateString (rho, A, T, E, I, L, s0, s1, k);

u1 = zeros(N1, 1);
u1Prev = zeros(N1, 1);
u1Next = zeros(N1, 1);

u2 = zeros(N1, 1);
u2Prev = zeros(N1, 1);
u2Next = zeros(N1, 1);

courantNo = c^2 * k^2 / h1^2 + 4 * kappa^2 * k^2 / h1^4
%% Bridge offset
b = 0.5;
barrier = -2e-6;

%% Collision Variables
cL1 = floor (N1 / 8); % bridge location
cL2 = floor (N2 / 2); % bridge location

%% Non-linear Spring Variables
K1 = 1000;
K3 = 0;
etaSpring = u1(cL1) - u2(cL2);
etaSpringPrev = u1(cL1) - u2(cL2);

%% Initial condition string because of bridge
% u1(3:cL) = (1:cL-2) / (cL - 2) * (b+2*eps);
% u1(cL:N-1) = (N-1-cL:-1:0) / (N-1-cL) * (b+2*eps);
% u1Prev = u1;

%% Excitation
width = 10;
loc = 1/5;
startIdx = floor(floor(loc * N1) - width / 2);
endIdx = floor(floor(loc * N1) + width / 2);
amp = 0.1;
u1(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
u1Prev = u1;

%% Initialise Energy Vectors
kinEnergy1 = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);
kinEnergy2 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);
connEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

rOCkinEnergy1 = zeros(lengthSound, 1);
rOCpotEnergy1 = zeros(lengthSound, 1);
rOCdamp0Energy = zeros(lengthSound, 1);
rOCdamp1Energy = zeros(lengthSound, 1);
rOCTotEnergy1 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCcolEnergy = zeros(lengthSound, 1);
rOCTotEnergy2 = zeros(lengthSound, 1);

rOCTotEnergy = zeros(lengthSound, 1);

Adiv = zeros(N1, 1);
v = zeros(N1, 1);

vec = 3:N1-2;
eVec = 2:N1-1; % energy vector
Jbr1 = zeros(N1,1);
Jbr1(cL1) = 1 / h1;

Jbr2 = zeros(N2,1);
Jbr2(cL2) = 1 / h2;

outputPos = floor(N1 / 6);

qPrev=0;
etaSpringPrev = 0;
etaSpring = 0;
for n = 2:lengthSound

    %% Update FDS
    
    etaSpring = u1Next(cL1) - u2Next(cL2);
    
    u1Next = B1 * u1 + C1 * u1Prev;
    u2Next = B2 * u2 + C2 * u2Prev;
   
%     Fphi = K1 + 2 * K3 * etaSpring^2;
%     if Fphi == 0
%         Falpha = 0;
%     else
%         Falpha = (2 * K1 / Fphi * etaSpring + etaSpringPrev + u1Next(cL) - u2Next)...
%             / (4 / Fphi + 1 / (h * (rho * A / k^2 + s0/k)) + rho * A * h * k^2 / M^2);
% u1NextTest(vec) = (rho * A / k^2 * (2 * u1(vec) - u1Prev(vec)) + T / h^2 * (u1(vec+1) - 2 * u1(vec) + u1(vec-1))...
%     - E * I / h^4 * (u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2))) * k^2 / (rho * A);
%     Falpha = (2 * etaSpring + etaSpringPrev + u1Next(cL1) - u2Next(cL2)) ...
%       / (4 / K1 + k^2 /(h1 * (rho * A)) + k^2 / (h2 * rho * A))
    varPhi = K1/4 + (K3 * etaSpring^2)/2;
%     Falpha = (-u1Next(cL1) + u2Next(cL2) * varPhi - (K1 * etaSpring) / 2 - etaSpringPrev * varPhi) / (2/h1 * varPhi + 1);
%  
%     Falpha = ((K1 / 2 * etaSpring) / varPhi + etaSpringPrev + u1Next(cL1) - u2Next(cL1))...
%         / (1 / varPhi + k^2 * rho * A / h1 + k^2 * rho * A / h2);
    Falpha = (u1Next(cL1) - u2Next(cL2) + 2 * etaSpring + etaSpringPrev)...
        / (4 / K1 + k^2 / (h1 * rho * A) + k^2 / (h2 * rho * A));
%     end
    u1Next = u1Next - (Jbr1 * Falpha) * (k^2 / (rho * A));
    u2Next = u2Next + (Jbr2 * Falpha) * (k^2 / (rho * A));

    %% Calculate energy of the system
    kinEnergy1(n) = rho * A / 2 * h1 * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h1 * sum((u1(3:N1) - u1(2:N1-1)) .* (u1Prev(3:N1) - u1Prev(2:N1-1)))...
        + E * I / 2 * 1/h1^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    kinEnergy2(n) = rho * A / 2 * h2 * sum((1/k * (u2 - u2Prev)).^2);
    potEnergy2(n) = T / 2 * 1/h2 * sum((u2(3:N2) - u2(2:N2-1)) .* (u2Prev(3:N2) - u2Prev(2:N2-1)))...
        + E * I / 2 * 1/h1^3 * sum((u2(eVec+1) - 2 * u2(eVec) + u2(eVec-1)) ...
        .* (u2Prev(eVec+1) - 2 * u2Prev(eVec) + u2Prev(eVec-1)));
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    connEnergy(n) = K1 / 2 * (1/2 * (etaSpring + etaSpringPrev))^2;
    
 
%     colEnergy(n) = 1/2 * psiPrev^2;
    totEnergy(n) = energy1(n) + energy2(n) + connEnergy(n);

%% Update states
%     psiPrevMS = psiMS; 
%     psiPrev = psiMB; 
    u1Prev = u1;
    u1 = u1Next;
%     u1PrevTest = u1Test;
%     u1Test = u1NextTest;
    u2Prev = u2;
    u2 = u2Next;
    
    etaSpringPrev = etaSpring;
    out(n) = u1Next(outputPos);
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(3,1,1);
        cla
        plot(u1Next, 'Linewidth', 2);
        hold on;
        plot(u2Next, 'Linewidth', 2);
        ylim([-0.1 0.1])
%         plot([cL-1, cL+1], [barrier barrier], 'Linewidth', 5);
%         ylim([min(u1Next), max(u1Next)])
        grid on; 
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
%         if ~strcmp(excitation, "bowed")
%             ylim([-amp, amp])
%         end
%         legend(["String", "Bridge", "Body"])
        subplot(3,1,2);
        cla
%         plot(energy1(10:n))
%         hold on;
%         plot(energy2(10:n))
%         plot(connEnergy(10:n))
% plot(energy1(10:n) + connEnergy(10:n) + energy2(10:n))
        hold on;
                 hold on;
        plot(connEnergy(10:n))
        plot(energy1(10:n))
%         hold on;
        plot(energy2(10:n))
        subplot(3,1,3)
        cla

%         plot(energy1(10:n))
        plot(totEnergy(10:n) / totEnergy(10) - 1)
hold on
% plot(energy2(10:n))
% plot(energy1(10:n))
% plot(u1Next - u1NextTest)
%         cla
%         plot(rOCenergy1(10:n));
%         hold on;
%         plot(rOCenergy2(10:n))
        drawnow
    end
end
if ~drawThings
    plot(out)
end