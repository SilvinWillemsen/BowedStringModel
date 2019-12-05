clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Excitation (cos or bowed)
exc = "bowed";
plotWeirdBehaviour = true;

%% Drawing Functions
drawThings = true;
drawSpeed = 10;
lengthSound = fs * 2;
drawStart = fs*0.4;
damping = true;
dampTest = false;
onlyString = false;

%% Bridge offset and location
offset = 1e-5;
bridgeLoc = 0.8;

%% String Variables
f0 = 100;
rhoS = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
% T = 10;  % Tension
T = c^2 * rhoS * A;
k = 1 / fs;         % Time step
ES = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappaS = sqrt (ES*I / (rhoS*A));   % Stiffness coefficient

% Damping coefficients
if damping
    s0S = 0.1;
    s1S = 0.05;
else
    s0S = 0;
    s1S = 0;
end

[BS, CS, NS, hS, Dxx, Dxxxx, s0S, s1S, bB, bC] = unscaledCreateStringNR (rhoS, A, T, ES, I, L, s0S, s1S, k);

u1 = zeros(NS, 1) + offset;
u1Prev = zeros(NS, 1) + offset;
u1Next = zeros(NS, 1) + offset;

% e = ones(NS-1,1);
% Dxx = (1/hS^2)*spdiags([e -2*e e],-1:1,NS-1,NS-1);
% 
courantNoS = c^2 * k^2 / hS^2 + 4 * kappaS^2 * k^2 / hS^4

%% Bowing terms
bP = floor (1/pi * NS);
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
u2Next = 0;

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
Nx = Nx - 1;
Ny = Ny - 1;

u3 = zeros(NP, 1);
horPos = 0.6;
vertPos = bridgeLoc;
horOutPos = 1 / 3;
vertOutPos = 1 / 7;

brP = floor(Ny * floor(horPos * Nx) + Ny * vertPos);
outputPosPlate = floor(Ny * floor(horOutPos * Nx) + Ny * vertOutPos);


if exc == "plate"
    u3(brP) = -offset;
end

u3Prev = zeros(NP, 1);
u3Next = zeros(NP, 1);

courantNoP = kappaP * k / hP^2

%% Collision Variables
cL = floor (NS * bridgeLoc); % bridge location
alpha1 = 1.3;
K1 = 5 * 10^10;

alpha2 = 1.3;
K2 = 5 * 10^10;


%% Excitation
amp = -offset*10;

if strcmp(exc, "cos")
    width = 10;
    loc = 1/4;
    startIdx = floor(floor(loc * NS) - width / 2);
    endIdx = floor(floor(loc * NS) + width / 2);
    u1Next(startIdx : endIdx) = u1(startIdx : endIdx) - amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
end
u1 = u1Next;
u1Prev = u1Next;

% amp = 10 * offset;
% if strcmp(exc, "cos")
%     width = 0.1 * NS;
%     loc = 1/4;
%     startIdx = floor(floor(loc * NS) - width / 2);
%     endIdx = floor(floor(loc * NS) + width / 2);
%     u1(startIdx : endIdx) = u1(startIdx : endIdx) - amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
% end
% u1Prev = u1;

%% Initialise
psi1Prev = 0;
psi2Prev = 0;

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

colEnergy1 = zeros(lengthSound, 1);
colEnergy2 = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

%% Initialise Rate-Of-Change Vectors (used to check damping)
rOCkinEnergy1 = zeros(lengthSound, 1);
rOCpotEnergy1 = zeros(lengthSound, 1);
rOCcolEnergy1 = zeros(lengthSound, 1);
rOCdamp0StringEnergy = zeros(lengthSound, 1);
rOCdamp1StringEnergy = zeros(lengthSound, 1);
rOCenergy1 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCcolEnergy2 = zeros(lengthSound, 1);
rOCenergy2 = zeros(lengthSound, 1);

rOCkinEnergy3 = zeros(lengthSound, 1);
rOCpotEnergy3 = zeros(lengthSound, 1);
rOCcolEnergy3 = zeros(lengthSound, 1);
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
for n = 1:lengthSound 
%     if n < rampVal
%         Fb = n * FbInit / rampVal;
%     else
        Fb = FbInit;
%     end
    eta1 = u2 - u1(cL);
    eta2 = u3(brP) - u2;
    %% Explicitly calculate g for Mass-Plate
%     if alpha1 == 1
%         g1 = 0;
%         if eta1 > 0
%             g1 = sqrt(K1 * (alpha1+1) / 2);
%         end
%     else
        g1 = sqrt(K1 * (alpha1+1) / 2) * abs(eta1)^((alpha1 - 1)/2);
%     end
    g1Save(n) = g1;
    
%     if g1 ~= 0
%         disp("wait")
%     end
    if alpha2 == 1
        g2 = 0;
        if eta2 > 0
            g2 = sqrt(K2 * (alpha2+1) / 2);
        end
    else
        g2 = sqrt(K2 * (alpha2+1) / 2) * subplus(eta2)^((alpha2 - 1)/2);
    end
    g2Save(n) = g2;
    
    %% Update FDSs without connection-force and collision terms 
    strVec = 3:NS-2;
    
    u1Next(strVec) = BS(strVec, :) * u1 + CS(strVec, :) * u1Prev;
%     u1Next = 2*u1-u1Prev+T*k^2 /(rhoS*A) * Dxx*u1;
    u2Next = (2 * u2 - u2Prev - k^2 * w1^2 * (u2 - offset) + R * k / 2 * u2Prev) / (1 + R * k / 2);
%     u2NextTest = ((2-w1^2*k^2) * u2 + (-1 + R * k / (2 * M)) * u2Prev) / (1 + R * k / (2 * M));
    u3Next = BP * u3 + CP * u3Prev;
    
    %% Bowing
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
    
    s0S = s0S / (rhoS * A);
    s1S = s1S / (rhoS * A);
    s0P = s0P / (rhoP * H);
    s1P = s1P / (rhoP * H);
    
    %% Matrix solution
    
    Adiv = [1 + g1^2*k^2 / (2*M*(2+R*k)) + g1^2*k^2 / (4*rhoS*A*hS*(1+s0S*k)), -g2^2*k^2 / (2*M*(2+R*k)); ...
            -g1^2*k^2 / (2*M*(2+R*k)), 1 + g2^2*k^2 / (2*M*(2+R*k)) + g2^2*k^2 / (4*hP^2*rhoP*H*(1+s0P*k)) ];
        
    v = [(u2Next - u2Prev - u1Next(cL) + u1Prev(cL)) / (2*k) - k*(psi1Prev*g1 - psi2Prev*g2) / (M*(2+R*k)) - psi1Prev*g1*k / (2*rhoS*A*hS*(1+s0S*k)); ...
        (u3Next(brP) - u3Prev(brP) - u2Next + u2Prev) / (2*k) + k*(psi1Prev*g1 - psi2Prev*g2) / (M*(2+R*k)) - psi2Prev*g2*k / (2*rhoP*H*hP^2*(1+s0P*k)) ];
    
    solut = Adiv \ v;
    psi1 = g1*k*solut(1) + psi1Prev;
    psi2 = g2*k*solut(2) + psi2Prev;
    
%     u1Next = u1Next + (- Jbow * Fb * BM * q * exp(-a*q^2)) ./ ((rhoS * A) / k^2 + s0S / k);
    u1Next = u1Next + k^2 / (rhoS * A * (1 + s0S*k)) * (Jbr * 0.5 * (psi1 + psi1Prev) * g1 - Jbow * Fb * BM * q * exp(-a*q^2));
    u2Next = u2Next - k^2 / (M * (1 + R*k / 2)) * 0.5 * (psi1 + psi1Prev) * g1 + k^2 / (M * (1 + R*k / 2)) * 0.5 * (psi2 + psi2Prev) * g2;
    u3Next = u3Next - JbrP * k^2 / (rhoP * H * (1 + s0P*k)) * 0.5 * (psi2 + psi2Prev) * g2;
%     Falpha = 0;
    s0S = s0S * (rhoS * A);
    s1S = s1S * (rhoS * A);
    s0P = s0P * (rhoP * H);
    s1P = s1P * (rhoP * H);
%%%%% undamped stuff
%     Adiv = [1 + g1^2*k^2 / (4*M) + g1^2*k^2 / (4*rhoS*A*hS), -g2^2*k^2 / (4*M); ...
%             -g1^2*k^2 / (4*M), 1 + g2^2*k^2 / (4*M) + g2^2*k^2 / (4*hP^2*rhoP*H) ];
%         
%     v = [(u2Next - u2Prev - u1Next(cL) + u1Prev(cL)) / (2*k) - k*(psi1Prev*g1 - psi2Prev*g2) / (2*M) - psi1Prev*g1*k / (2*rhoS*A*hS); ...
%         (u3Next(brP) - u3Prev(brP) - u2Next + u2Prev) / (2*k) + k*(psi1Prev*g1 - psi2Prev*g2) / (2*M) - psi2Prev*g2*k / (2*rhoP*H*hP^2) ];
%     
%     solut = Adiv \ v;
%     
%     psi1 = g1*k*solut(1) + psi1Prev;
%     psi2 = g2*k*solut(2) + psi2Prev;
%     
%     u1Next = u1Next + Jbr * k^2 / (rhoS * A) * 0.5 * (psi1 + psi1Prev) * g1;
%     u2Next = u2Next - k^2 / M * 0.5 * (psi1 + psi1Prev) * g1 + k^2 / (M * (1 + R*k / 2)) * 0.5 * (psi2 + psi2Prev) * g2;
%     u3Next = u3Next - JbrP * k^2 / (rhoP * H) * 0.5 * (psi2 + psi2Prev) * g2;
    if drawThings && drawStart == 0
%         disp("Difference between etaNext and u3Next(brP) - u2Next (should be 0)")
%         etaNext - (u3Next(brP) - u2Next)
    elseif mod(n,floor(lengthSound / 100)) == 0
        prog = prog + 1;
        disp("Progress: " + prog + "%")
    end
    
%     etaSpringNext = u1Next(cL) - u2Next;
    
    %% Update Psi
%     psi1 = psi1Prev + 0.5 * g1 .* (eta1Next - eta1Prev);
%     psi2 = psi2Prev + 0.5 * g2 .* (eta2Next - eta2Prev);
    
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
    colEnergy1(n) = psi1Prev^2 / 2;
    colEnergy2(n) = psi2Prev^2 / 2;
    
    % Total Energy
    totEnergy(n) = energy1(n) + energy2(n) + energy3(n) + colEnergy1(n) + colEnergy2(n);
    
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
    rOCdampMassEnergy(n) = -M * R / (4 * k^2) * (u2Next - u2Prev) * (u2Next - u2Prev);
    rOCenergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n) - rOCdampMassEnergy(n);
    
    % Plate
    rOCkinEnergy3(n) = hP^2 * rhoP * H / (2*k^3) * sum((u3Next - 2 * u3 + u3Prev) .* (u3Next - u3Prev));
    rOCpotEnergy3(n) = -(D*hP^2) / (2*k) * sum((DD * u3) .* (u3Next - u3Prev));
    rOCdamp0PlateEnergy(n) = -2 * s0P * hP^2 / (4 * k^2) * sum((u3Next - u3Prev) .* (u3Next - u3Prev));
    rOCdamp1PlateEnergy(n) = 2 * s1P * hP^2 / (2 * k^2 * hP^2) * sum((Denergy * u3 - Denergy * u3Prev) .* (u3Next - u3Prev));
    
    rOCenergy3(n) = rOCkinEnergy3(n) - rOCpotEnergy3(n) - rOCdamp0PlateEnergy(n) - rOCdamp1PlateEnergy(n);
   
    % Collision (string part)
    rOCcolEnergy1(n) = 1 / (4 * k) * g1 * (psi1 + psi1Prev) * (u1Next(cL) - u1Prev(cL));
    % Collision (mass part)
    rOCcolEnergy2(n) = -1 / (4 * k) * g1 * (psi1 + psi1Prev) * (u2Next - u2Prev) + 1 / (4 * k) * g2 * (psi2 + psi2Prev) * (u2Next - u2Prev);
    % Collision (plate part)
    rOCcolEnergy3(n) = -1 / (4 * k) * g2 * (psi2 + psi2Prev) * (u3Next(brP) - u3Prev(brP));
        
    %Total Energy
    rOCTotEnergy(n) = rOCenergy1(n) + rOCenergy2(n) + rOCenergy3(n) - rOCcolEnergy1(n) - rOCcolEnergy2(n) - rOCcolEnergy3(n); %including damping so should be 0
    
    %% Update states
    psi1Prev = psi1; 
    psi2Prev = psi2; 
    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    u3Prev = u3;
    u3 = u3Next;
%     
%     eta1Prev = eta1;
%     eta1 = eta1Next;
%     
%     eta2Prev = eta2;
%     eta2 = eta2Next;
    
    qPrev = q;
    
    out(n) = u1Next(outputPos);
    out2(n) = u3Next(outputPosPlate);
    out3(n) = u2Next;
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        if ~plotWeirdBehaviour
            % Draw States of... 
            subplot(4,1,1);
            cla
            hold on;
            %...string
            plot(u1Next, 'Linewidth', 2);    
            %...mass
            scatter(cL, u2Next, 400, '.');    
            %...barrier
            plot([cL-1, cL+1], [u3Next(brP) u3Next(brP)], 'Linewidth', 5);  

            % Extra functions
    %             ylim([-amp / 2, amp / 2]); % Set y-limit to the amplitude of the raised cosine
            grid on; 
            set(gca, 'Linewidth', 2, 'Fontsize', 16)
            title("State of the system")
            legend(["String", "Bridge", "Body"])

            %...plate
            subplot(4,1,2)
            imagesc(reshape(u3Next, [Ny,Nx])')
            xticks([])
            yticks([])
            title("Body (plate)")
            set(gca,'Fontsize', 16)
            subplot(4,1,3);
            if s0S == 0 && s1S == 0 && s0P == 0 && s1P == 0
                % Draw Normalised energy
                plot(totEnergy(10:n) / totEnergy(10) - 1);
    %             hold off;
    %             plot((kinEnergy1(1:n) + potEnergy1(1:n)) - (kinEnergy1(1) + potEnergy1(1)))
    %             hold on;
    %             plot(colEnergy1(1:n))
    %             plot(energy3(10:n) / energy3(10) - 1);
                title("Normalised Energy")
            elseif exc == "cos"
                plot(totEnergy(10:n))
                title("Total Energy (should decrease)")
            else
                range = 100;
                if (n > range)
                    plot(g2Save(n-range:n));
                end
                title("$g_2^n$", 'interpreter', 'latex')
    %             title("Total Energy (should be variable due to bow)");
            end
            subplot(4,1,4)
            % Draw rate of change of the energy
            cla
            hold on;
    %             plot(rOCdamp0PlateEnergy(10:n))
    %         plot(rOCbowStringEnergy(1:n))
            plot(rOCTotEnergy(1:n))
    %         plot(g1Save(1:n))
    %         plot(rOCbowStringEnergy(1:n))
            title("Rate of change of Energy minus damping (should be 0 within machine precision)")

    %         subplot(3,1,3)
    %         plot(g1Save(1:n));
    %         title("$g^n$", 'interpreter', 'latex')
    %         set(gca,'Fontsize', 16)
        else
            subplot(5,1,1)
            plot(energy1(1:n));
            subplot(5,1,2)
            plot(energy2(1:n));
            subplot(5,1,3)
            plot(energy3(1:n));
            subplot(5,1,4)
            plot(colEnergy1(1:n));
            subplot(5,1,5)
            plot(colEnergy2(1:n));
        end
        drawnow;
    end
end
posOut1 = out - min(out);
totOut1 = (posOut1/max(abs(posOut1)) - 0.5) * 2;
posOut2 = out2 - min(out2);
totOut2 = (posOut2/max(abs(posOut2)) - 0.5) * 2;
posOut3 = out3 - min(out3);
totOut3 = (posOut3/max(abs(posOut3)) - 0.5) * 2;
totOut4 = (g1Save/max(abs(g1Save)) - 0.5) * 2;
totOut = totOut1/5 + totOut2 + totOut3;
% soundsc(totOut1 /5 + totOut2 + totOut3,fs)

if ~drawThings
    if onlyString
        plot(out);
    else
        plot(totOut);
    end
end
plot(g2Save)