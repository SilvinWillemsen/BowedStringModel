clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Excitation (cos or bowed)
exc = "bowed";

%% Drawing Functions
drawThings = false;
drawSpeed = 10000;
lengthSound = fs*3;
drawStart = 0;
damping = true;
dampTest = false;
onlyString = false;

%% Bridge offset and location
offset = 1e-5;
bridgeLoc = 1 / 20;

%% String Variables
f0 = 100;    
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
if damping
    s0S = 0.1;
    s1S = 0.005;
else
    s0S = 0;
    s1S = 0;
end

[BS, CS, NS, hS, Dxx, Dxxxx, s0S, s1S, bB, bC] = unscaledCreateStringNR (rhoS, A, T, ES, I, L, s0S, s1S, k);

u1 = zeros(NS, 1) + offset;
u1Prev = zeros(NS, 1) + offset;
u1Next = zeros(NS, 1) + offset;

courantNoS = c^2 * k^2 / hS^2 + 4 * kappaS^2 * k^2 / hS^4

%% Bowing terms
bP = floor (0.5 * NS);
a = 100;
BM = sqrt(2 * a) * exp(1/2);

if strcmp(exc, "bowed")
    FbInit = 1;
else
    FbInit = 0;
end

Vb = -0.2;
qPrev = -Vb;    
tol = 1e-4;

%% Mass Variables
f1 = 1000;    % fundamental frequency [Hz] (< 1 / (k * pi) (< 14,037 Hz))
w1 = 2 * pi * f1;   % angular frequency (< 2 / k (< 88,200 rad/s))
M = 0.001;

if damping
    R = 0.1;
else
    R = 0.0;
end

u2 = offset;
u2Prev = u2;
u2Next = u2;

%% Plate Variables
Lx = 0.2;
Ly = 1.5;
rhoP = 7850;
EP = 2e11;
H = 0.001;

if damping
    s0P = 5;
    s1P = 0.01;
else
    s0P = 0;
    s1P = 0;
end

[BP, CP, NP, Nx, Ny, hP, kappaP, Denergy, D, DD, s0P, s1P] = unscaledCreatePlate (Lx, Ly, rhoP, EP, H, s0P, s1P, k);

u3 = zeros(NP, 1);
horPos = 0.6;
vertPos = bridgeLoc;
horOutPos = 1 / 3;
vertOutPos = 1 / 7;

brP = floor((Ny - 1) * floor(horPos * (Nx - 1)) + (Ny - 1) * vertPos);
outputPosPlate = floor((Ny - 1) * floor(horOutPos * (Nx - 1)) + (Ny - 1) * vertOutPos);

% u3(brP) = offset;
u3Prev = zeros(NP, 1);
u3Next = zeros(NP, 1);

courantNoP = kappaP * k / hP^2

%% Collision Variables
cL = floor (NS * bridgeLoc); % bridge location
alpha = 1.3;
K = 5 * 10^10;

%% Non-linear Spring Variables
if ~onlyString
    K1 = 1000000;
    K3 = 100;
    if damping
        sx = 1;
    else
        sx = 0;
    end
else
    K1 = 0;
    K3 = 0;
    sx = 0;
end

etaSpring = u1(cL) - u2;
etaSpringPrev = u1(cL) - u2;

%% Excitation
amp = 100 * offset;
if strcmp(exc, "cos")
    width = 10;
    loc = 1/4;
    startIdx = floor(floor(loc * NS) - width / 2);
    endIdx = floor(floor(loc * NS) + width / 2);
    u1(startIdx : endIdx) = u1(startIdx : endIdx) - amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
end
u1Prev = u1;

%% Initialise
etaPrev = u3(brP) - u2;
psiPrev = 0;
eta = u3(brP) - u2Prev;
etaNext = u3(brP) - u2;

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

%% Initialise Rate-Of-Change Vectors (used to check damping)
rOCkinEnergy1 = zeros(lengthSound, 1);
rOCpotEnergy1 = zeros(lengthSound, 1);
rOCdamp0StringEnergy = zeros(lengthSound, 1);
rOCdamp1StringEnergy = zeros(lengthSound, 1);
rOCenergy1 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCcolEnergy1 = zeros(lengthSound, 1);
rOCenergy2 = zeros(lengthSound, 1);

rOCkinEnergy3 = zeros(lengthSound, 1);
rOCpotEnergy3 = zeros(lengthSound, 1);
rOCcolEnergy2 = zeros(lengthSound, 1);
rOCenergy3 = zeros(lengthSound, 1);

rOCconnEnergy = zeros(lengthSound, 1);
rOCconnDampEnergy = zeros(lengthSound, 1);
rOCTotEnergy = zeros(lengthSound, 1);

vec = 3:NS-2;
eVec = 2:NS-1; % energy vector
Jbr = zeros(NS,1);
Jbr(cL) = 1 / hS;

JbrP = zeros(NP, 1);
JbrP(brP) = 1 / hP^2;
IbrP = zeros(NP, 1);
IbrP(brP) = 1;

Ibow = zeros(NS,1);
Ibow(bP) = 1;
Jbow = Ibow / hS;

outputPos = floor(NS / 6);
prog = 0;

rampVal = 10000;
for n = 2:lengthSound 
%     if n < rampVal
%         Fb = n * FbInit / rampVal;
%     else
        Fb = FbInit;
%     end
    
    %% Explicitly calculate g for Mass-Plate
    if alpha == 1
        g = 0;
        if eta > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end
    gSave(n) = g;
    
    %% Update FDSs without connection-force and collision terms 
    strVec = 3:NS-2;
    
    u1Next(strVec) = BS(strVec, :) * u1 + CS(strVec, :) * u1Prev;
    u2Next = (M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * (u2 - offset) + R / (2*k) * u2Prev) / (M / k^2 + R / (2*k));
    u3Next = BP * u3 + CP * u3Prev;
    
    if exc == "bowed"
        b = 2/k * Vb + 2 * s0S * Vb + Ibow' * bB * u1 + Ibow' * bC * u1Prev;
        eps = 1;
        i = 0;

        while eps>tol && i < 100
            q=qPrev-(1/(rhoS * A * hS) * Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0S*qPrev+b)/...
             (1/(rhoS * A * hS)*Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0S);
            eps = abs(q-qPrev);
            qPrev = q;
            i = i + 1;
        end
    else
        q = 0;
    end
%     iSave(n) = i;
%     qSave(n) = q;
    
    
%     u1Next(strVec) = ((rhoS * A / k^2 + s0S / k) * (BS(strVec, :) * u1 + CS(strVec, :) * u1Prev)) ./ (rhoS * A / k^2 + s0S / k);
    
%     %% Calculate connection forces
    phiMinus = K1 / 4 + K3 * etaSpring^2 / 2 - sx / k;
    phiPlus = K1 / 4 + K3 * etaSpring^2 / 2 + sx / k;
    phiMinTest = k * (K1 + 2 * K3 * etaSpring^2) - 4 * sx;
    phiPlusTest = k * (K1 + 2 * K3 * etaSpring^2) + 4 * sx;
    varPsi = 1 / (hS * ((rhoS * A) / k^2 + s0S / k)) ...
        + (1/(M / k^2 + R / (2*k) + g^2/4)) + 1 / phiPlus;
%     varPsiTest = ((4 * k * hS * rhoS * A + k^2 * phiPlusTest) * (4 * M + g^2 * k^2) + (4 * k^2 * hS * rhoS * A * phiPlusTest)) ... 
%         / (hS * rhoS * A * phiPlusTest * (4 * M + g^2 * k^2));
    varPsiTest = phiPlusTest * k^2 * (4 * M + 2 * R * k + g^2 * k^2 + 4 * hS * (rhoS * A + s0S * k)) + 4 * k * (hS * (rhoS * A + s0S * k)) * (4 * M + 2 * R * k + g^2 * k^2);
    if phiPlus == 0
%         Falpha = 0;
        FalphaTick = 0;
    else
        % One division
        FalphaTick = (phiPlusTest * hS * (rhoS * A + s0S * k) * (4 * M + 2 * R * k + g^2 * k^2)...
            * ((4 * M + 2 * R * k + g^2 * k^2) * (2 * k * K1 * etaSpring + phiMinTest * etaSpringPrev + phiPlusTest * u1Next(cL)) ...
            - phiPlusTest * ((4 * M + 2 * R * k) * u2Next - g^2 * k^2 * etaPrev + 4 * k^2 * psiPrev * g)))...
            / (varPsiTest * phiPlusTest * (4 * M + 2 * R * k + g^2 * k^2));
    end
    
    plateTerm = (phiPlusTest * hS * (rhoS * A + s0S * k) * (4 * M + 2 * R * k + g^2 * k^2) * g^2 * k^2)...
        / (varPsiTest * (4 * M + 2 * R * k + g^2 * k^2));
    
    Adiv = [rhoS * A / k^2 + s0S / k,                 0,                            -plateTerm / hS;...
                     0,               M / k^2 + g^2 / 4 + R / (2*k),            plateTerm - g^2 / 4; ...
                     0,                       -g^2 / (4 * hP^2),      rhoP * H / k^2 + s0P / k + g^2 / (4*hP^2)];
    
    v = [... String
         rhoS * A / k^2 * (2 * u1(cL) - u1Prev(cL)) + T / hS^2 * (u1(cL+1) - 2 * u1(cL) + u1(cL-1)) ...
         - ES * I / hS^4 * (u1(cL+2) - 4 * u1(cL+1) + 6 * u1(cL) - 4 * u1(cL-1) + u1(cL-2)) - FalphaTick / hS...
         + s0S / k * u1Prev(cL)...
         + 2 * s1S / (hS^2 * k) * (u1(cL+1) - 2 * u1(cL) + u1(cL-1) - u1Prev(cL+1) + 2 * u1Prev(cL) - u1Prev(cL-1)); ...
         ... Mass
         M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * (u2 - offset) + R / (2*k) * u2Prev - g^2 / 4 * etaPrev + psiPrev * g + FalphaTick; ...
         ... Plate
         rhoP * H / k^2 * (2 * u3(brP) - u3Prev(brP)) - D * DD(brP, :) * u3...
         + s0P / k * u3Prev(brP)...
         + 2 * s1P / (hP^2 * k) * (Denergy(brP, :) * u3 - Denergy(brP, :) * u3Prev)...
         + (g^2 / 4 * etaPrev - psiPrev * g) / hP^2];
    
    solut = Adiv \ v;
    u2Next = solut(2);
    Falpha = FalphaTick - plateTerm * solut(3);
    etaNext = solut(3) - u2Next;
%     Falpha = 0;
    u1Next = (((rhoS * A) / k^2 + s0S / k) * u1Next - (Jbr * Falpha) - Jbow * Fb * BM * q * exp(-a*q^2)) ./ ((rhoS * A) / k^2 + s0S / k);
    u3Next = u3Next - JbrP * (g^2 / 4 * (etaNext - etaPrev) + psiPrev * g) / ((rhoP * H) / k^2 + s0P / k);
    
    if drawThings && drawStart == 0
%         disp("Difference between etaNext and u3Next(brP) - u2Next (should be 0)")
%         etaNext - (u3Next(brP) - u2Next)
    elseif mod(n,floor(lengthSound / 100)) == 0
        prog = prog + 1;
        disp("Progress: " + prog + "%")
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
    
    % String
    rOCkinEnergy1(n) = hS * rhoS * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = hS * T / (2*k*hS^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
         - hS * ES * I / (2 * k * hS^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0StringEnergy(n) = -2 * s0S * hS / (4 * k^2) * sum((u1Next - u1Prev).*(u1Next - u1Prev));
    rOCdamp1StringEnergy(n) = 2 * hS * s1S / (2 * k^2 * hS^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1) - u1Prev(eVec+1) + 2 * u1Prev(eVec) - u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    rOCbowStringEnergy(n) =  -Fb * (u1Next(bP) - u1Prev(bP)) / (2 * k) * sqrt(2*a) * q * exp(-a*q^2 + 1/2);
    %-sum(Jbow * hS * Fb * sqrt(2*a) .* ((u1Next - u1Prev) / (2 * k) - Vb) * exp(-a * q^2 + 1/2) .* (u1Next - u1Prev) / (2*k));
    rOCenergy1(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0StringEnergy(n) - rOCdamp1StringEnergy(n) - rOCbowStringEnergy(n);
    
    % Mass
    rOCkinEnergy2(n) = M / (2*k^3) * (u2Next - u2Prev) * (u2Next - 2 * u2 + u2Prev);
    rOCpotEnergy2(n) = -M * w1^2 / (2*k) * (u2Next - u2Prev) * (u2 - offset);
    rOCdampMassEnergy(n) = -R / (4 * k^2) * (u2Next - u2Prev) * (u2Next - u2Prev);
    rOCenergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n) - rOCdampMassEnergy(n);
    
    % Plate
    rOCkinEnergy3(n) = hP^2 * rhoP * H / (2*k^3) * sum((u3Next - 2 * u3 + u3Prev) .* (u3Next - u3Prev));
    rOCpotEnergy3(n) = -(D*hP^2) / (2*k) * sum((DD * u3) .* (u3Next - u3Prev));
    rOCdamp0PlateEnergy(n) = -2 * s0P * hP^2 / (4 * k^2) * sum((u3Next - u3Prev) .* (u3Next - u3Prev));
    rOCdamp1PlateEnergy(n) = 2 * s1P * hP^2 / (2 * k^2 * hP^2) * sum((Denergy * u3 - Denergy * u3Prev) .* (u3Next - u3Prev));
    
    rOCenergy3(n) = rOCkinEnergy3(n) - rOCpotEnergy3(n) - rOCdamp0PlateEnergy(n) - rOCdamp1PlateEnergy(n);

    % Collision (mass part)
    rOCcolEnergy1(n) = 1 / (4 * k) * g * (psi + psiPrev) * (u2Next - u2Prev);
    % Collision (plate part)
    rOCcolEnergy2(n) = -1 / (4 * k) * g * (psi + psiPrev) * (u3Next(brP) - u3Prev(brP));
    
    % Connection
    rOCconnEnergy(n) = (K1 / 4 * (etaSpringNext + 2 * etaSpring + etaSpringPrev) ...
        + K3 / 2 * etaSpring^2 * (etaSpringNext + etaSpringPrev))...
        * 1/(2*k) * (etaSpringNext - etaSpringPrev);
    rOCconnDampEnergy(n) = 2 * sx / (4 * k^2) * (etaSpringNext - etaSpringPrev)^2;
    
    %Total Energy
    rOCTotEnergy(n) = rOCenergy1(n) + rOCenergy2(n) + rOCenergy3(n) + rOCconnEnergy(n) + rOCconnDampEnergy(n) - rOCcolEnergy1(n) - rOCcolEnergy2(n); %including damping so should be 0
    
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
    
    qPrev = q;
    
    out(n) = u1Next(outputPos);
%     out2(n) = sum(u3Next);
    out2(n) = u3Next(outputPosPlate);
    out3(n) = u2Next;
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        
        % Draw States of... 
        
        subplot(3,1,1);
        cla
        hold on;
        %...string
        plot(u1Next, 'Linewidth', 2);    
        %...mass
        scatter(cL, u2Next, 400, '.');    
        %...barrier
        plot([cL-1, cL+1], [u3Next(brP) u3Next(brP)], 'Linewidth', 5);  
            
            % Extra functions
            ylim([-amp / 2, amp / 2]); % Set y-limit to the amplitude of the raised cosine
            grid on; 
            set(gca, 'Linewidth', 2, 'Fontsize', 16)
            title("State of the system")
            legend(["String", "Bridge", "Body"])
   
        %...plate
        subplot(3,1,2)
        imagesc(reshape(u3Next, [Ny-1,Nx-1])')
        xticks([])
        yticks([])
        title("Body (plate)")
        set(gca,'Fontsize', 16)
%         subplot(4,1,3);
%         if s0S == 0 && s1S == 0 && s0P == 0 && s1P == 0 && sx == 0
%             % Draw Normalised energy
%             plot(totEnergy(10:n) / totEnergy(10) - 1);
% %             plot(energy3(10:n) / energy3(10) - 1);
%             title("Normalised Energy")
%         elseif exc == "cos"
%             plot(totEnergy(10:n));
%             title("Total Energy (should decrease)")
%         else
%             plot(totEnergy(10:n));
%             title("Total Energy (should be variable due to bow)");
%         end
%         subplot(4,1,4)
%         % Draw rate of change of the energy
%         cla
%         hold on;
% %             plot(rOCdamp0PlateEnergy(10:n))
% %         plot(rOCbowStringEnergy(1:n))
%         plot(rOCTotEnergy(1:n))
% %         plot(rOCbowStringEnergy(1:n))
%         title("Rate of change of Energy minus damping (should be 0 within machine precision)")
%         
        subplot(3,1,3)
        plot(gSave(1:n));
        title("$g^n$", 'interpreter', 'latex')
        set(gca,'Fontsize', 16)
        drawnow;
    end
end
posOut1 = out - min(out);
totOut1 = (posOut1/max(abs(posOut1)) - 0.5) * 2;
posOut2 = out2 - min(out2);
totOut2 = (posOut2/max(abs(posOut2)) - 0.5) * 2;
posOut3 = out3 - min(out3);
totOut3 = (posOut3/max(abs(posOut3)) - 0.5) * 2;
totOut4 = (gSave/max(abs(gSave)) - 0.5) * 2;
totOut = totOut1/5 + totOut2 + totOut3;
% soundsc(totOut1 /5 + totOut2 + totOut3,fs)

if ~drawThings
    if onlyString
        plot(out);
    else
        plot(totOut);
    end
end
figure;
plot(gSave)