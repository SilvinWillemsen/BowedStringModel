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
s0 = 0.1;
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
K1 = 0;
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

vec = 3:N-2;
eVec = 2:N-1; % energy vector
Jbr = zeros(N,1);
Jbr(cL) = 1 / h;

outputPos = floor(N / 6);

for n = 2:lengthSound 
    %% Calculate g for Mass-Barrier
    if alpha == 1
        g = 0;
        if eta > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end
    
    %% Update FDSs
    u1Next = B * u1 + C * u1Prev;
    u2Next = M/k^2 * (2*u2 - u2Prev) - M * w1^2 * u2 + (barrier - etaPrev) * g^2 / 4 + psiPrev * g...
        / (M / k^2 + g^2 / 4);
    
    %% Add connection forces
    if K1 == 0 && K3 == 0
        Falpha = 0;
    else
        varPhi = K1 / 4 + K3 * etaSpring^2 / 2;
        Falpha = (u1Next(cL) - u2Next + K1 * etaSpring / (2 * varPhi) + etaSpringPrev) ...
             / (k^2 / (h * (rho * A)) + (1/(M / k^2 + g^2/4)) + 1 / varPhi);
    end
    u1Next = u1Next - (Jbr * Falpha) * (k^2 / (rho * A));
    u2Next = u2Next + Falpha / (M/k^2 + g^2/4);
    
    %% Update eta
    etaNext = barrier - u2Next;
 
    %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
    
    %% Calculate energy of the system
    % Energy String
    kinEnergy1(n) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h * sum((u1(3:N) - u1(2:N-1)) .* (u1Prev(3:N) - u1Prev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    % Energy Mass
    kinEnergy2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2(n) = M / 2 * w1^2 * (u2 * u2Prev);
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    % Connection and Collision energies
    connEnergy(n) = K1 / 2 * (1/2 * (etaSpring + etaSpringPrev))^2 + ...
        K3 / 4 * (etaSpring^2 * etaSpringPrev^2);
    colEnergy(n) = psiPrev^2 / 2;
    
    % Total Energy
    totEnergy(n) = energy1(n) + energy2(n) + connEnergy(n) + colEnergy(n);
    
    %% Calculate rate-of-changes in energies (inner product of delta tdot with scheme) 
    rOCkinEnergy1(n) = h * rho * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = h * T / (2*k*h^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
         - h * E * I / (2 * k * h^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0Energy(n) = -rho * A * 2 * s0 * h / (4 * k^2) * sum((u1Next - u1Prev).^2);
    rOCdamp1Energy(n) = rho * A * 2 * h * s1 / (2 * k^2 * h^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1) - u1Prev(eVec+1) + 2 * u1Prev(eVec) - u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    rOCcolEnergy(n) = h / (4 * k) * sum(g .* (psi + psiPrev) .* (u1Next - u1Prev));
    rOCTotEnergy(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0Energy(n) - rOCdamp1Energy(n) - rOCcolEnergy(n); %including damping so should be 0
    
    %% Update states
    psiPrev = psi; 
    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    etaSpringPrev = etaSpring;
    etaSpring = u1Next(cL) - u2Next;
    
    etaPrev = eta;
    eta = etaNext;
    
    out(n) = u1Next(outputPos);
    
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        % Draw States of 
        subplot(2,1,1);
        cla
        hold on;
        plot(u1Next, 'Linewidth', 2);                           % String
        scatter(cL, u2Next, 400, '.');                          % Mass
        plot([cL-1, cL+1], [barrier barrier], 'Linewidth', 5);  % Barrier
        
        % Set y-limit to the amplitude of the raised cosine
        ylim([-amp, amp])
        
        % Extra functions
        grid on; 
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        title("State of the system")
        legend(["String", "Mass", "Barrier"])
        
        subplot(2,1,2);
        if s0 == 0 && s1 == 0
            % Draw Normalised energy
            plot(totEnergy(10:n) / totEnergy(10) - 1);
            title("Normalised Energy")
        else 
            % Draw rate of change of the energy
            plot(rOCTotEnergy(10:n))
            title("Rate of change of Energy minus damping")
        end
        drawnow
    end
end
if ~drawThings
    plot(out)
end