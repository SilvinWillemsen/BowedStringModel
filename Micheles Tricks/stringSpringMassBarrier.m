clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = true;
drawSpeed = 1;
lengthSound = 2 * fs;
drawStart = 1;

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

[B, C, N, h, Dxx, Dxxxx] = unscaledCreateString (rho, A, T, E, I, L, s0, s1, k);

u1 = zeros(N, 1);
u1Prev = zeros(N, 1);
u1Next = zeros(N, 1);

courantNo = c^2 * k^2 / h^2 + 4 * kappa^2 * k^2 / h^4
%% Bridge offset
b = 0.5;
barrier = -2e-6;

%% Mass Variables
f1 = 0;            % fundamental frequency [Hz]
w1 = 2 * pi * f1;   % angular frequency
M = 0.001;
u2 = 0.1;
u2Prev = 0.1;
% u2Next = -0.1;

%% Collision Variables
cL = floor (N / 8); % bridge location

alpha = 1.3;
K = 0 * 10^8;

%% Non-linear Spring Variables
K1 = 1000;
K3 = 0;
etaSpring = u1(cL) - u2;
etaSpringPrev = u1(cL) - u2;

%% Initial condition string because of bridge
% u1(3:cL) = (1:cL-2) / (cL - 2) * (b+2*eps);
% u1(cL:N-1) = (N-1-cL:-1:0) / (N-1-cL) * (b+2*eps);
% u1Prev = u1;

%% Excitation
% width = 10;
% loc = 1/5;
% startIdx = floor(floor(loc * N) - width / 2);
% endIdx = floor(floor(loc * N) + width / 2);
% amp = 0.1;
% u1(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
% u1Prev = u1;

%% Initialise
etaPrev = barrier - u2;
psiPrev = 0;
eta = barrier - u2Prev;
etaNext = barrier - u2;

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

Adiv = zeros(N, 1);
v = zeros(N, 1);

vec = 3:N-2;
eVec = 2:N-1; % energy vector
Jbr = zeros(N,1);
Jbr(cL) = 1 / h;


outputPos = floor(N / 6);

qPrev=0;
MassRatio = rho * A * h / M;
u1NextTest = zeros(N,1);
for n = 2:lengthSound
    %% Calculate energy of the system
    kinEnergy1(n) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h * sum((u1(3:N) - u1(2:N-1)) .* (u1Prev(3:N) - u1Prev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    kinEnergy2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2(n) = M / 2 * w1^2 * (u2 * u2Prev);
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    connEnergy(n) = K1 / 2 * (1/2 * (etaSpring + etaSpringPrev))^2;
    
 
%     colEnergy(n) = 1/2 * psiPrev^2;
    totEnergy(n) = energy1(n) + energy2(n) + connEnergy(n);
    
    
    %% Calculate g for Mass Barrier
    if alpha == 1
        g = 0;
        if eta > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end
    
%     %% Matrix Form
%     Amat = [rho*A/k^2 + gMS.^2/(4*h) + (rho * A * s0)/k,...
%             -gMS.^2/(4*h);...
%             -gMS.^2/4,...
%             M/k^2 + gMS.^2/4 + g.^2/4];
%     answ = [rho*A/k^2 * (2 * u1(cL) - u1Prev(cL)) + T / h^2 * (u1(cL+1) - 2 * u1(cL) + u1(cL-1)) ...
%             - E * I / h^4 * (u1(cL+2) - 4 * u1(cL+1) + 6 * u1(cL) - 4 * u1(cL-1) + u1(cL-2))...
%             + rho * A * s0 / k * u1Prev(cL)...
%             + 2 * rho * A * s1 / (k * h^2) * (u1(cL+1) - 2 * u1(cL) + u1(cL-1) - u1Prev(cL+1) + 2 * u1Prev(cL) - u1Prev(cL-1))...
%             + (-gMS.^2 / 4 * etaPrevMS + psiPrevMS .* gMS) / h;...
%             M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2 + gMS.^2/4 * etaPrevMS - psiPrevMS .* gMS...
%             + g.^2/4 * u2Prev + psiPrev * g];
%     solut = Amat\answ;
%     
%     etaNext = barrier - solut(2);
%     u2Next = solut(2);
%     
    %% Update FDS
    u1Next = B * u1 + C * u1Prev;
    u2Next = (M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2) * k^2 / M;
    etaSpringPrev = etaSpring;
    etaSpring = u1Next(cL) - u2Next;
%     Fphi = K1 + 2 * K3 * etaSpring^2;
%     if Fphi == 0
%         Falpha = 0;
%     else
%         Falpha = (2 * K1 / Fphi * etaSpring + etaSpringPrev + u1Next(cL) - u2Next)...
%             / (4 / Fphi + 1 / (h * (rho * A / k^2 + s0/k)) + rho * A * h * k^2 / M^2);
% u1NextTest(vec) = (rho * A / k^2 * (2 * u1(vec) - u1Prev(vec)) + T / h^2 * (u1(vec+1) - 2 * u1(vec) + u1(vec-1))...
%     - E * I / h^4 * (u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2))) * k^2 / (rho * A);
      Falpha = (u1Next(cL) - u2Next + 2 * etaSpring + etaSpringPrev) ...
          / (k^2 /(h * (rho * A)) + k^2 * MassRatio / M + 4 / K1);
%     end
    u1Next = u1Next - (Jbr * Falpha) * (k^2 / (rho * A));
    u2Next = u2Next + (Falpha * MassRatio) * (k^2 / M);
%     u1Next(vec) = (rho * A / k^2 * (2 * u1(vec) - u1Prev(vec)) + T / h^2 * (u1(vec+1) - 2 * u1(vec) + u1(vec-1)) ...
%          - E * I / h^4 * (u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2))...
%          + rho * A * s0 / k * u1Prev(vec)...
%          + 2 * rho * A * s1 / (k * h^2) * (u1(vec+1) - 2 * u1(vec) + u1(vec-1) - u1Prev(vec+1) + 2 * u1Prev(vec) - u1Prev(vec-1))...
%          + Jbr(vec) * (gMS^2/4 * (etaNextMS - etaPrevMS) + psiPrevMS * gMS)...
%          + exciteSignal) ./ (rho * A / k^2 + rho * A * s0/k + Jbow(vec) * Btot / (2*k));
%      u1NextTest(vec) = (2 * u1Test(vec) - u1PrevTest(vec) + lambdaSq * (u1Test(vec+1) - 2 * u1Test(vec) + u1Test(vec-1))...
%          - muSq *  (u1Test(vec+2) - 4 * u1Test(vec+1) + 6 * u1Test(vec) - 4 * u1Test(vec-1) + u1Test(vec-2))...
%          + s0 * k * u1PrevTest(vec)...
%          + 2 * s1 * k / h^2 * (u1Test(vec+1) - 2 * u1Test(vec) + u1Test(vec-1) - u1PrevTest(vec+1) + 2 * u1PrevTest(vec) - u1PrevTest(vec-1))) / (1 + s0 * k);
%          
%     if drawThings
%         etaNextMS - (u2Next - u1Next(cL));
%         etaNext - (b - u2Next);
%     end
%     %% Update Psi
%     psiMS = psiPrevMS + 0.5 * gMS .* (etaNextMS - etaPrevMS);
%     psiMB = psiPrev + 0.5 * g .* (etaNext - etaPrev);
%     
% %     %% Calculate rate-of-changes in energy for checking damping (inner product of delta t-dot with scheme)     
    rOCkinEnergy1(n) = h * rho * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = h * T / (2*k*h^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
         - h * E * I / (2 * k * h^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0Energy(n) = -2 * rho * A * s0 * h / (4 * k^2) * sum((u1Next - u1Prev).^2);
    rOCdamp1Energy(n) = 2 * h * rho * A * s1 / (2 * k^2 * h^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1) - u1Prev(eVec+1) + 2 * u1Prev(eVec) - u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    rOCenergy1(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0Energy(n) - rOCdamp1Energy(n);
%     rOCcolEnergyMS1(n) = gMS / (4 * k) * (psiMS + psiPrevMS) * (u1Next(cL) - u1Prev(cL));
%     rOCcolEnergyMS2(n) = -gMS / (4 * k) * (psiMS + psiPrevMS) * (u2Next - u2Prev);
%     
    rOCkinEnergy2(n) = M / (2 * k^3) * (u2Next - 2 * u2 + u2Prev) * (u2Next - u2Prev);
    rOCpotEnergy2(n) = -M * w1^2 / (2*k) * (u2Next - u2Prev) * u2;
    rOCenergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n);
%     
%     rOCcolEnergy(n) = g / (4 * k) * (psiMB + psiPrev) * (u2Next - u2Prev);
    rOCTotEnergy(n) = rOCenergy1(n) + rOCenergy2(n); % - rOCcolEnergyMS1(n) - rOCcolEnergyMS2(n) - rOCcolEnergy(n); %including damping so should be 0
    %% Update states
%     psiPrevMS = psiMS; 
%     psiPrev = psiMB; 
    u1Prev = u1;
    u1 = u1Next;
%     u1PrevTest = u1Test;
%     u1Test = u1NextTest;
    u2Prev = u2;
    u2 = u2Next;
    out(n) = u1Next(outputPos);
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(3,1,1);
        cla
        plot(u1Next, 'Linewidth', 2);
        hold on;
        scatter(cL, u2Next, 400, '.');
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
        plot(connEnergy(10:n))
        hold on;
        plot(energy1(10:n))
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