clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = true;
drawSpeed = 100;
lengthSound = 2 * fs;
drawStart = 0;

%% Offset
offset = 1e-3;

%% String Variables
f0 = 0;    
rhoS = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rhoS * A;  % Tension
k = 1 / fs;         % Time step
ES = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappaS = sqrt (ES*I / (rhoS*A));   % Stiffness coefficient

% Damping coefficients
s0S = 0.0;
s1S = 0.000;

[BS, CS, NS, hS, Dxx, Dxxxx] = unscaledCreateString (rhoS, A, T, ES, I, L, s0S, s1S, k);

u1 = zeros(NS, 1) + offset;
u1Prev = zeros(NS, 1) + offset;
u1Next = zeros(NS, 1) + offset;

courantNoS = c^2 * k^2 / hS^2 + 4 * kappaS^2 * k^2 / hS^4

%% Mass Variables
f1 = 0;            % fundamental frequency [Hz] (doesn't work with offset)
w1 = 2 * pi * f1;   % angular frequency
M = 0.001;
u2 = offset;
u2Prev = u2;
u2Next = u2;

%% Plate Variables
Lx = 0.5;
Ly = 2;
rhoP = 7850;
EP = 2e11;
H = 0.005;
s0P = 0;
s1P = 0;

[BP, CP, NP, Nx, Ny, hP, kappaP, Denergy, D, DD] = unscaledCreatePlate (Lx, Ly, rhoS, EP, H, s0P, s1P, k);
plateLoc = 0; % body offset

u3 = zeros(NP, 1);
brP = floor(NP / 2 + Nx / 2);
% u3(brP) = 1;
u3Prev = zeros(NP, 1);
u3Next = zeros(NP, 1);

courantNoP = kappaP * k / hP^2

%% Collision Variables
cL = floor (NS / 8); % bridge location
alpha = 1.3;
K = 5 * 10^10;

%% Non-linear Spring Variables
K1 = 100000;
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
startIdx = floor(floor(loc * NS) - width / 2);
endIdx = floor(floor(loc * NS) + width / 2);
amp = 20*offset;
u1(startIdx : endIdx) = u1(startIdx : endIdx) - amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
u1Prev = u1;

%% Initialise
etaPrev = plateLoc - u2;
psiPrev = 0;
eta = plateLoc - u2Prev;
etaNext = plateLoc - u2;

%% Initialise Energy Vectors
kinEnergy1 = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);

kinEnergy2 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);

kinEnergy3 = zeros(lengthSound, 1);
potEnergy3 = zeros(lengthSound, 1);
energy3 = zeros(lengthSound, 1);

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
rOCcolEnergy1 = zeros(lengthSound, 1);
rOCTotEnergy2 = zeros(lengthSound, 1);

rOCTotEnergy = zeros(lengthSound, 1);

vec = 3:NS-2;
eVec = 2:NS-1; % energy vector
Jbr = zeros(NS,1);
Jbr(cL) = 1 / hS;

JbrP = zeros(NP, 1);
JbrP(brP) = 1 / hP^2;
IbrP = zeros(NP, 1);
IbrP(brP) = 1;

outputPos = floor(NS / 6);

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
    strVec = 3:NS-2;
    %% Update FDSs
    u1Next(strVec) = BS(strVec, :) * u1 + CS(strVec, :) * u1Prev;
    u2Next = (M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2) * k^2 / M;
    u3Next = BP * u3 + CP * u3Prev;
    
    %% Add connection forces
    varPhi = K1 / 4 + K3 * etaSpring^2 / 2;
    varPsi = k^2 / (hS * (rhoS * A)) + (1/(M / k^2 + g^2/4)) + 1 / varPhi;
    if varPhi == 0
%         Falpha = 0;
        FalphaTick = 0;
    else
%         Falpha = (u1Next(cL) - u2Next + K1 * etaSpring / (2 * varPhi) + etaSpringPrev) / varPsi;
        FalphaTick = (K1 * etaSpring / (2 * varPhi) + etaSpringPrev + u1Next(cL)...
            -(M/k^2 * u2Next - g^2/4 * etaPrev + psiPrev * g) / (M / k^2 + g^2 / 4)) / varPsi;
    end
    plateTerm = (g^2/4)/(varPsi * (M / k^2 + g^2 / 4));
    
    Adiv = [rhoS * A / k^2 + s0S / k,         0,                  -plateTerm;...
                   0,                 M / k^2 + g^2 / 4,      plateTerm - g^2 / 4; ...
                   0,                 -g^2 / (4 * hP^2),   rhoP * H / k^2 + g^2 / (4*hP^2)];
    
    v = [rhoS * A / k^2 * (2 * u1(cL) - u1Prev(cL)) + T / hS^2 * (u1(cL+1) - 2 * u1(cL) + u1(cL-1)) ...
         - ES * I / hS^4 * (u1(cL+2) - 4 * u1(cL+1) + 6 * u1(cL) - 4 * u1(cL-1) + u1(cL-2)) - FalphaTick / hS...
         + s0S / k * u1Prev(cL)...
         + 2 * s1S / (hS^2 * k) * (u1(cL+1) - 2 * u1(cL) + u1(cL-1) - u1Prev(cL+1) + 2 * u1Prev(cL) - u1Prev(cL-1)); ...
         M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * (u2 - offset) - g^2 / 4 * etaPrev + psiPrev * g + FalphaTick; ...
         rhoP * H / k^2 * (2 * u3(brP) - u3Prev(brP)) - D * DD(brP, :) * u3 + (g^2 / 4 * etaPrev - psiPrev * g) / hP^2];
    
    solut = Adiv \ v;
    u2Next = solut(2);
    Falpha = FalphaTick - plateTerm * solut(3);
    etaNext = solut(3) - u2Next;
    u1Next = u1Next - (Jbr * Falpha) / ((rhoS * A) / k^2 + s0S / k);
    u3Next = u3Next - JbrP * (g^2 / 4 * (etaNext - etaPrev) + psiPrev * g) * (k^2 / (rhoP * H));
    
%     Adiv = [M/k^2 + g^2/4, -g^2/4; ...
%             -g^2/(4 * hP^2), rhoP * H / k^2 + g^2 / (4 * hP^2)];
%     v = [M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2 + Falpha - g^2 / 4 * etaPrev + psiPrev * g; ...
%          rhoP * H / k^2 * (2 * u3(brP) - u3Prev(brP)) - D * DD(brP, :) * u3 + (g^2 / 4 * etaPrev - psiPrev * g) / hP^2];
% 
%     solut = Adiv\v;
%     etaNext = solut(2) - solut(1);
%     u2Next = solut(1); 
%     
%     u3Next = u3Next - JbrP * (g^2 / 4 * (etaNext - etaPrev) + psiPrev * g) * (k^2 / (rhoP * H));
    
    if drawThings
        etaNext - (u3Next(brP) - u2Next)
    end
    
    etaSpringNext = u1Next(cL) - u2Next;
    
    %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
    
    %% Calculate energy of the system
    % Energy String
    kinEnergy1(n) = rhoS * A / 2 * hS * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/hS * sum((u1(3:NS) - u1(2:NS-1)) .* (u1Prev(3:NS) - u1Prev(2:NS-1)))...
        + ES * I / 2 * 1/hS^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    % Energy Mass
    kinEnergy2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2(n) = M / 2 * w1^2 * ((u2 - offset) * (u2Prev - offset));
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    
    % Energy Plate
    kinEnergy3(n) = ((rhoP * H) / 2) * hP^2 * sum(sum(1/k^2 * (u3 - u3Prev).^2));
    potEnergy3(n) = D / (2 * hP^2) * sum((Denergy * u3) .* (Denergy * u3Prev));
    energy3(n) = kinEnergy3(n) + potEnergy3(n);
    
    % Connection and Collision energies
    connEnergy(n) = K1 / 2 * (1/2 * (etaSpring + etaSpringPrev))^2 + ...
        K3 / 4 * (etaSpring^2 * etaSpringPrev^2);
    colEnergy(n) = psiPrev^2 / 2;
    
    % Total Energy
    totEnergy(n) = energy1(n) + energy2(n) + energy3(n) + connEnergy(n) + colEnergy(n);
    
    %% Calculate rate-of-changes in energies (inner product of delta tdot with scheme) 
    rOCkinEnergy1(n) = hS * rhoS * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = hS * T / (2*k*hS^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
         - hS * ES * I / (2 * k * hS^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0Energy(n) = -rhoS * A * 2 * s0S * hS / (4 * k^2) * sum((u1Next - u1Prev).^2);
    rOCdamp1Energy(n) = rhoS * A * 2 * hS * s1S / (2 * k^2 * hS^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1) - u1Prev(eVec+1) + 2 * u1Prev(eVec) - u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    rOCenergy1(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0Energy(n) - rOCdamp1Energy(n);
    
    rOCkinEnergy2(n) = M / (2*k^3) * (u2Next - u2Prev) * (u2Next - 2 * u2 + u2Prev);
    rOCpotEnergy2(n) = -M * w1^2 / (2*k) * (u2Next - u2Prev) * u2;
    rOCenergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n);

    rOCcolEnergy1(n) = 1 / (4 * k) * sum(g * (psi + psiPrev) .* (u2Next - u2Prev));
    rOCcolEnergy2(n) = -hP^2 / (4 * k) * sum(g * (psi + psiPrev) .* (u3Next(brP) - u3Prev(brP)));
    
    rOCconnEnergy(n) = (K1 / 4 * (etaSpringNext + 2 * etaSpring + etaSpringPrev) ...
        + K3 / 2 * etaSpring^2 * (etaSpringNext + etaSpringPrev))...
        * 1/(2*k) * (etaSpringNext - etaSpringPrev);
    
    rOCTotEnergy(n) = rOCenergy1(n) + rOCenergy2(n) + rOCconnEnergy(n); %including damping so should be 0
    
    %% Update states
    psiPrev = psi; 
    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    u3Prev = u3;
    u3 = u3Next;
    
    etaSpringPrev = etaSpring;
    etaSpring = u1Next(cL) - u2Next;
    
    etaPrev = eta;
    eta = etaNext;
    
    out(n) = u1Next(outputPos);
    
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        % Draw States of 
        subplot(3,2,1);
        cla
        hold on;
        plot(u1Next, 'Linewidth', 2);                           % String
        scatter(cL, u2Next, 400, '.');                          % Mass
%         plot([cL-1, cL+1], [plateLoc plateLoc], 'Linewidth', 5);  % Barrier
        subplot(3,2,[2, 4, 6])
        imagesc(reshape(u3Next, [Ny-1,Nx-1]))
        % Set y-limit to the amplitude of the raised cosine
%         ylim([-amp, amp])
        
        % Extra functions
%         grid on; 
%         set(gca, 'Linewidth', 2, 'Fontsize', 16)
%         title("State of the system")
%         legend(["String", "Mass", "Barrier"])
        
        subplot(3,2,3);
%         if s0S == 0 && s1S == 0
%             % Draw Normalised energy
            plot(totEnergy(10:n) / totEnergy(10) - 1);
% %             plot(energy3(10:n) / energy3(10) - 1);
%             title("Normalised Energy")
%         else 
            % Draw rate of change of the energy
%             plot(rOCTotEnergy(10:n))
%             title("Rate of change of Energy minus damping")
%         end
        subplot(3,2,5)
        cla;
        hold on
%         plot(energy1(10:n))
%         plot(energy2(10:n))
%         plot(energy3(10:n))
        plot(colEnergy(10:n))
        drawnow
    end
end
if ~drawThings
    plot(out)
end