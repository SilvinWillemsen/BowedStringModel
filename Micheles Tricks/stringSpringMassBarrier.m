clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = true;
drawSpeed = 1000;
lengthSound = 2 * fs;
drawStart = 0;

%% String Variables
f0 = 100;    
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rho * A;  % Tension
k = 1 / fs;         % Time` step
E = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappa = sqrt (E*I / (rho*A));   % Stiffness coefficient

% Damping coefficients
s0 = 0.0;
s1 = 0.000;

[B, C, N, h, Dxx, Dxxxx] = unscaledCreateString (rho, A, T, E, I, L, s0, s1, k);

u1 = zeros(N, 1);
u1Prev = zeros(N, 1);
u1Next = zeros(N, 1);

courantNo = c^2 * k^2 / h^2 + 4 * kappa^2 * k^2 / h^4
%% Bridge offset
b = 0.0;
barrier = -2e-6;

%% Mass Variables
f1 = 0;            % fundamental frequency [Hz]
w1 = 2 * pi * f1;   % angular frequency
M = 0.1;
u2 = 0.0;
u2Prev = u2;
u2Next = u2;

%% Collision Variables
cL = floor (N / 8); % bridge location
alpha = 1.3;
K = 5 * 10^10;

%% Non-linear Spring Variables
K1 = 10000000;
K3 = 0;
etaSpring = u1(cL) - u2;
etaSpringPrev = u1(cL) - u2;

%% Initial condition string because of bridge
% u1(3:cL) = (1:cL-2) / (cL - 2) * (b+2*eps);
% u1(cL:N-1) = (N-1-cL:-1:0) / (N-1-cL) * (b+2*eps);
% u1Prev = u1;

%% Excitation
width = 10;
loc = 1/5;
startIdx = floor(floor(loc * N) - width / 2);
endIdx = floor(floor(loc * N) + width / 2);
amp = 0.1;
u1(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
u1Prev = u1;

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
colEnergy = zeros(lengthSound, 1);
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

% Adiv = zeros(N, 1);
% v = zeros(N, 1);

vec = 3:N-2;
eVec = 2:N-1; % energy vector
Jbr = zeros(N,1);
Jbr(cL) = 1 / h;

outputPos = floor(N / 6);

qPrev=0;
u1NextTest = zeros(N,1);
for n = 2:lengthSound
    etaSpring = u1Next(cL) - u2Next;
%     energy2(n) = M * (1/k * (u2 - u2Prev))^2 / 2 + M * w1^2*u2*u2Prev / 2 + psiPrev^2 / 2;
    etaPrev = eta;
    eta = barrier - u2Next;
    %% Calculate g for Mass Barrier
    if alpha == 1
        g = 0;
        if eta > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end
    gSave(n) = g;
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
    Adiv = M / k^2 + g^2 / 4;
    v = M/k^2 * (2*u2 - u2Prev) - M * w1^2 * u2 + (barrier - etaPrev) * g^2 / 4 + psiPrev * g;
    u2Next = v / Adiv;
    
    %% Update FDS
    u1Next = B * u1 + C * u1Prev;
    
%     Fphi = K1 + 2 * K3 * etaSpring^2;
    if K1 == 0 && K3 == 0
        Falpha = 0;
    else
%         Falpha = (2 * K1 / Fphi * etaSpring + etaSpringPrev + u1Next(cL) - u2Next)...
%             / (4 / Fphi + 1 / (h * (rho * A / k^2 + s0/k)) + rho * A * h * k^2 / M^2);
% u1NextTest(vec) = (rho * A / k^2 * (2 * u1(vec) - u1Prev(vec)) + T / h^2 * (u1(vec+1) - 2 * u1(vec) + u1(vec-1))...
%     - E * I / h^4 * (u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2))) * k^2 / (rho * A);
      
        varPhi = K1 / 4 + K3 * etaSpring^2 / 2;
%     Falpha = (u1Next(cL) - u2Next + 2 * etaSpring + etaSpringPrev) ...
%               / (k^2 / (h * (rho * A)) + (1/(M / k^2 + g^2/4)) + 4 / K1);
        Falpha = (u1Next(cL) - u2Next + K1 * etaSpring / (2 * varPhi) + etaSpringPrev) ...
            / (k^2 / (h * (rho * A)) + (1/(M / k^2 + g^2/4)) + 1 / varPhi);
    
    end
    u1Next = u1Next - (Jbr * Falpha) * (k^2 / (rho * A));
    u2Next = u2Next + Falpha / (M/k^2 + g^2/4);
    etaNext = barrier - u2Next;
    
    %% Calculate energy of the system
    kinEnergy1(n) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h * sum((u1(3:N) - u1(2:N-1)) .* (u1Prev(3:N) - u1Prev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    kinEnergy2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2(n) = M / 2 * w1^2 * (u2 * u2Prev);
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    connEnergy(n) = K1 / 2 * (1/2 * (etaSpring + etaSpringPrev))^2 + ...
        K3 / 4 * (etaSpring^2 * etaSpringPrev^2);
    colEnergy(n) = psiPrev^2 / 2;
    
    totEnergy(n) = energy1(n) + energy2(n) + connEnergy(n) + colEnergy(n);
%     totEnergy(n) = energy2(n);
    %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
     
    %% Update states
    psiPrev = psi; 
    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    etaSpringPrev = etaSpring;
    out(n) = u1Next(outputPos);
    
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(2,1,1);
        cla
        plot(u1Next, 'Linewidth', 2);
        hold on;
        scatter(cL, u2Next, 400, '.');
        plot([cL-1, cL+1], [barrier barrier], 'Linewidth', 5);
        ylim([-amp, amp])
        grid on; 
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        subplot(2,1,2);
%          plot(totEnergy(10:n))
        plot(totEnergy(10:n) / totEnergy(10) - 1);
        drawnow
    end
end
if ~drawThings
    plot(out)
end