clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Excitation (cos or bowed)
exc = "cos";
plotWeirdBehaviour = false;

%% Drawing Functions
drawThings = true;
drawSpeed = 5;
lengthSound = fs * 2;
drawStart = 0;
damping = false;
dampTest = false;
onlyString = false;

%% Bridge offset and location
offset = 1e-5;
bridgeLoc = 0.8;

%% String Variables
f0 = 50;
rhoS = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
% T = 10;  % Tension
T = c^2 * rhoS * A;
k = 1 / fs;         % Time step
ES = 0;           % Young's modulus
Iner = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappaS = sqrt (ES*Iner / (rhoS*A));   % Stiffness coefficient

% Damping coefficients
if damping
    s0S = 0.1;
    s1S = 1;
else
    s0S = 0;
    s1S = 0;
end

[BS, CS, NS, hS, Dxx, Dxxxx, s0S, s1S, bB, bC] = unscaledCreateStringNR (rhoS, A, T, ES, Iner, L, s0S, s1S, k);

u1 = zeros(NS, 1) + offset;
u1Prev = zeros(NS, 1) + offset;
u1Next = zeros(NS, 1) + offset;

w1 = zeros(NS, 1) + offset;
w1Prev = zeros(NS, 1) + offset;
w1Next = zeros(NS, 1) + offset;

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
quPrev = -Vb;    
qwPrev = -Vb;
tol = 1e-4;

%% Mass Variables
f1 = 100;    % fundamental frequency [Hz] (< 1 / (k * pi) (< 14,037 Hz))
omega0 = 2 * pi * f1;   % angular frequency (< 2 / k (< 88,200 rad/s))
M = 0.001;

if damping
    R = 0.1;
else
    R = 0.0;
end

u2 = offset;
u2Prev = u2 * 1.00001;
% u2Prev = -1e-10;
u2Next = 0;

w2 = offset;
w2Prev = w2 * 1.00001;
% w2Prev = -1e-10;
w2Next = 0;

%% Plate Variables
Lx = 1.5;
Ly = 0.4;
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
w3 = zeros(NP, 1);

horPos = 0.5;
vertPos = 0.5;
horOutPos = 1 / 3;
vertOutPos = 1 / 7;

brP = floor(Ny * floor(horPos * Nx) + Ny * vertPos);
outputPosPlate = floor(Ny * floor(horOutPos * Nx) + Ny * vertOutPos);

if exc == "plate"
    u3(brP) = -offset;
    w3(brP) = -offset;
end

u3Prev = u3;
u3Next = zeros(NP, 1);
w3Prev = w3;
w3Next = zeros(NP, 1);

courantNoP = kappaP * k / hP^2

%% Collision Variables
cL = floor (NS * bridgeLoc) + 1; % bridge location
alpha1 = 1.3;
K1 = 5 * 10^10;

alpha2 = 1.3;
K2 = 5 * 10^10;

%% Excitation
% amp = -offset*10000;
amp = -1e-4;
if strcmp(exc, "cos")
    width = floor(NS / 5);
    loc = 28/40;
    startIdx = floor(floor(loc * NS) - width / 2) + 1;
    endIdx = floor(floor(loc * NS) + width / 2) + 1;
    u1Next(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
    w1Next(startIdx : endIdx) = w1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
end
u1 = u1Next;
u1Prev = u1Next;

w1 = w1Next;
w1Prev = w1Next;

%% Initialise
psi1PrevU = 0;
psi2PrevU = 0;
psi1PrevW = 0;
psi2PrevW = 0;

%% Initialise Energy Vectors U
kinEnergy1u = zeros(lengthSound, 1);
potEnergy1u = zeros(lengthSound, 1);
energy1u = zeros(lengthSound, 1);

kinEnergy2u = zeros(lengthSound, 1);
potEnergy2u = zeros(lengthSound, 1);
energy2u = zeros(lengthSound, 1);

kinEnergy3u = zeros(lengthSound, 1);
potEnergy3u = zeros(lengthSound, 1);
energy3u = zeros(lengthSound, 1);

colEnergy1u = zeros(lengthSound, 1);
colEnergy2u = zeros(lengthSound, 1);
totEnergyu = zeros(lengthSound, 1);

%% Initialise Energy Vectors W
kinEnergy1w = zeros(lengthSound, 1);
potEnergy1w = zeros(lengthSound, 1);
energy1w = zeros(lengthSound, 1);

kinEnergy2w = zeros(lengthSound, 1);
potEnergy2w = zeros(lengthSound, 1);
energy2w = zeros(lengthSound, 1);

kinEnergy3w = zeros(lengthSound, 1);
potEnergy3w = zeros(lengthSound, 1);
energy3w = zeros(lengthSound, 1);

colEnergy1w = zeros(lengthSound, 1);
colEnergy2w = zeros(lengthSound, 1);
totEnergyw = zeros(lengthSound, 1);

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
DeltaDelta = Denergy * Denergy;

etawNext1 = w2 - w1(cL);
etawNext2 = w3(brP) - w2;

colTol = 1e-23;
etaSolut = [0, 0];
for n = 1:lengthSound 
    
    Fb = FbInit;
        
    eta1u = u2 - u1(cL);
    eta2u = u3(brP) - u2;
    %% Explicitly calculate g for String-Mass

    g1 = sign(eta1u) * sqrt(K1 * (alpha1+1) / 2) * abs(eta1u)^((alpha1 - 1)/2);
    g1Save(n) = g1;
    
    g2 = 0;
    if alpha2 == 1
        if eta2u > 0
            g2 = sqrt(K2 * (alpha2+1) / 2);
        end
    else
        g2 = sqrt(K2 * (alpha2+1) / 2) * subplus(eta2u)^((alpha2 - 1)/2);
    end
    g2Save(n) = g2;
    
    %% Update FDSs without connection-force and collision terms 
    strVec = 3:NS-2;
    
    u1Next(strVec) = BS(strVec, :) * u1 + CS(strVec, :) * u1Prev;
    u2Next = (2 * u2 - u2Prev - k^2 * omega0^2 * (u2 - offset) + R * k / 2 * u2Prev) / (1 + R * k / 2);
    u3Next = BP * u3 + CP * u3Prev;
    
    %% Bowing
    if exc == "bowed"
        b = 2/k * Vb + 2 * s0S * Vb + Ibow' * bB * u1 + Ibow' * bC * u1Prev;
        eps = 1;
        i = 0;

        while eps>tol && i < 100
            qu=quPrev-(1/(rhoS * A * hS) * Fb*BM*quPrev*exp(-a*quPrev^2)+2*quPrev/k+2*s0S*quPrev+b)/...
             (1/(rhoS * A * hS)*Fb*BM*(1-2*a*quPrev^2)*exp(-a*quPrev^2)+2/k+2*s0S);
            eps = abs(qu-quPrev);
            quPrev = qu;
            i = i + 1;
        end
    else
        qu = 0;
    end
    s0S = s0S / (rhoS * A);
    s1S = s1S / (rhoS * A);
    s0P = s0P / (rhoP * H);
    s1P = s1P / (rhoP * H);
    
    u1Next = u1Next - k^2 / (rhoS * A * (1 + s0S*k)) * (Jbow * Fb * BM * qu * exp(-a*qu^2));

    %% Matrix solution
    
    Adiv = [1 + g1^2*k^2 / (2*M*(2+R*k)) + g1^2*k^2 / (4*rhoS*A*hS*(1+s0S*k)), -g2^2*k^2 / (2*M*(2+R*k)); ...
            -g1^2*k^2 / (2*M*(2+R*k)), 1 + g2^2*k^2 / (2*M*(2+R*k)) + g2^2*k^2 / (4*hP^2*rhoP*H*(1+s0P*k)) ];
        
    v = [(u2Next - u2Prev - u1Next(cL) + u1Prev(cL)) / (2*k) - k*(psi1PrevU*g1 - psi2PrevU*g2) / (M*(2+R*k)) - psi1PrevU*g1*k / (2*rhoS*A*hS*(1+s0S*k)); ...
        (u3Next(brP) - u3Prev(brP) - u2Next + u2Prev) / (2*k) + k*(psi1PrevU*g1 - psi2PrevU*g2) / (M*(2+R*k)) - psi2PrevU*g2*k / (2*rhoP*H*hP^2*(1+s0P*k)) ];

%     a11 = 1.0 + g1^2 * k^2 / (2*M*(2+R*k)) + g1^2*k^2 / (4*rhoS*A*hS*(1+s0S*k));
%     a12 = -g2^2*k^2 / (2*M*(2+R*k));
%     a21 = -g1^2*k^2 / (2*M*(2+R*k));
%     a22 = 1.0 + g2^2*k^2 / (2*M*(2+R*k)) + g2^2*k^2 / (4*hP^2*rhoP*H*(1+s0P*k));
%  
%     oOdet = 1.0 / (a11 * a22 - a12 * a21);
%     
%     solut1Test = (v(1) * a22 - v(2) * a12) * oOdet;
%     solut2Test = (-v(1) * a21 + v(2) * a11) * oOdet;
    
    solut = Adiv \ v;
    
    psi1U = g1*k*solut(1) + psi1PrevU;
    psi2U = g2*k*solut(2) + psi2PrevU;
    
%     u1Next = u1Next + (- Jbow * Fb * BM * q * exp(-a*q^2)) ./ ((rhoS * A) / k^2 + s0S / k);
    u1Next = u1Next + k^2 / (rhoS * A * (1 + s0S*k)) * (Jbr * 0.5 * (psi1U + psi1PrevU) * g1);
    u2Next = u2Next - (k^2 / (M * (1 + R*k / 2))) * 0.5 * (psi1U + psi1PrevU) * g1 + (k^2 / (M * (1 + R*k / 2))) * 0.5 * (psi2U + psi2PrevU) * g2;
    u3Next = u3Next - JbrP * k^2 / (rhoP * H * (1 + s0P*k)) * 0.5 * (psi2U + psi2PrevU) * g2;
%     Falpha = 0;
    s0S = s0S * (rhoS * A);
    s1S = s1S * (rhoS * A);
    s0P = s0P * (rhoP * H);
    s1P = s1P * (rhoP * H);

    if drawThings && drawStart == 0
%         disp("Difference between etaNext and u3Next(brP) - u2Next (should be 0)")
%         etaNext - (u3Next(brP) - u2Next)
    elseif mod(n,floor(lengthSound / 100)) == 0
        prog = prog + 1;
        disp("Progress: " + prog + "%")
    end
    
    %% Newton Raphson
    etaw1 = w2 - w1(cL);
    etawPrev1 = w2Prev - w1Prev(cL);
    
    etaw2 = w3(brP) - w2;
    etawPrev2 = w3Prev(brP) - w2Prev;
    
    b1 = -2 * etaw1 + etawPrev1 + k^2 * omega0^2 * (w2 - offset) ...
    + T * k^2 / (rhoS * A * hS^2) * (w1(cL+1) - 2 * w1(cL) + w1(cL-1))...
    - ES * Iner * k^2 / (rhoS * A * hS^4) * (w1(cL+2) - 4 * w1(cL+1) + 6 * w1(cL) - 4 * w1(cL-1) + w1(cL-2));
    
    b2 = -2 * w3(brP) + w3Prev(brP) + 2 * w2 - w2Prev + D * k^2 / (rhoP * H * hP^4) * DeltaDelta(brP, :) * w3 - k^2 * omega0^2 * (w2 - offset);
    
    eps = 1;
    ii = 0;

    while eps > colTol && ii < 100
        
        % terms to use for the quotient rule
        f1 = abs(etawNext1)^(alpha1+1) - abs(etawPrev1)^(alpha1+1);
        df1 = (alpha1+1) * abs(etawNext1)^alpha1 * sign(etawNext1);
        gg1 = etawNext1 - etawPrev1;
        dgg1 = 1;
        
        f2 = subplus(etawNext2)^(alpha2+1) - subplus(etawPrev2)^(alpha2+1);
        df2 = (alpha2+1) * subplus(etawNext2)^alpha2 * 0.5 * (1 + sign(etawNext2));
        gg2 = etawNext2 - etawPrev2;
        dgg2 = 1;
        
        coeff1 = (k^2 / (rhoS * A * hS) + k^2 / M) * K1 / (alpha1+1);
        coeff2 = (k^2 / (rhoP * H * hP^2) + k^2 / M) * K2 / (alpha2+1);
        
        F = etawNext1 + coeff1 * f1 / gg1 - k^2 / M * K2 / (alpha2+1) * f2 / gg2 + b1;
        G = etawNext2 + coeff2 * f2 / gg2 - k^2 / M * K1 / (alpha1+1) * f1 / gg1 + b2;
        
        dF1 = 1 + coeff1 * (df1 * gg1 - f1 * dgg1) / gg1^2;
        dF2 = -k^2 / M * K2 / (alpha2+1) * (df2 * gg2 - f2 * dgg2 / gg2^2);
        
        dG1 = -k^2 / M * K1 / (alpha1+1) * (df1 * gg1 - f1 * dgg1 / gg1^2);
        dG2 = 1 + coeff2 * (df2 * gg2 - f2 * dgg2) / gg2^2;
        
        Jac = [dF1, dF2; dG1, dG2];
        prevEtaSolut = etaSolut;
        etaSolut = [etawNext1; etawNext2] - Jac \ [F; G];
        
        etawNextNew1 = etaSolut(1);
        etawNextNew2 = etaSolut(2);
        
        eps = norm(etaSolut - prevEtaSolut);
        etawNext1 = etawNextNew1;
        etawNext2 = etawNextNew2;
        
        ii = ii + 1;
    end
    disp(ii)
    
    force1(n) = K1 / (alpha1 + 1) * (abs(etawNext1)^(alpha1+1) - abs(etawPrev1)^(alpha1+1)) / (etawNext1 - etawPrev1);
    force2(n) = K2 / (alpha2 + 1) * (subplus(etawNext2)^(alpha2+1) - subplus(etawPrev2)^(alpha2+1)) / (etawNext2 - etawPrev2);
    
    w1Next(strVec) = BS(strVec, :) * w1 + CS(strVec, :) * w1Prev + Jbr(strVec) * k^2 / (rhoS * A) * force1(n);
    w2Next = (2 * w2 - w2Prev - k^2 * omega0^2 * (w2 - offset) + R * k / 2 * w2Prev) / (1 + R * k / 2) + k^2 / M * (-force1(n) + force2(n));
    w3Next = BP * w3 + CP * w3Prev - JbrP * k^2 / (rhoP * H) * force2(n);

    %% Calculate energy of the system
    % Energy String
    kinEnergy1u(n) = rhoS * A / 2 * hS * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1u(n) = T / 2 * 1/hS * sum((u1(3:NS) - u1(2:NS-1)) .* (u1Prev(3:NS) - u1Prev(2:NS-1)))...
        + ES * Iner / 2 * 1/hS^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1u(n) = kinEnergy1u(n) + potEnergy1u(n);
    
    % Energy Mass
    kinEnergy2u(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2u(n) = M / 2 * omega0^2 * ((u2 - offset) * (u2Prev - offset));
    energy2u(n) = kinEnergy2u(n) + potEnergy2u(n);
    
    % Energy Plate
    kinEnergy3u(n) = ((rhoP * H) / 2) * hP^2 * sum(sum(1/k^2 * (u3 - u3Prev).^2));
    potEnergy3u(n) = D / (2 * hP^2) * sum((Denergy * u3) .* (Denergy * u3Prev));
    energy3u(n) = kinEnergy3u(n) + potEnergy3u(n);
    
    % Connection and Collision energies
    colEnergy1u(n) = psi1PrevU^2 / 2;
    colEnergy2u(n) = psi2PrevU^2 / 2;
    
    % Total Energy
    totEnergyu(n) = energy1u(n) + energy2u(n) + energy3u(n) + colEnergy1u(n) + colEnergy2u(n);
    
    %% Calculate energy of the system
    % Energy String
    kinEnergy1w(n) = rhoS * A / 2 * hS * sum((1/k * (w1 - w1Prev)).^2);
    potEnergy1w(n) = T / 2 * 1/hS * sum((w1(3:NS) - w1(2:NS-1)) .* (w1Prev(3:NS) - w1Prev(2:NS-1)))...
        + ES * Iner / 2 * 1/hS^3 * sum((w1(eVec+1) - 2 * w1(eVec) + w1(eVec-1)) ...
        .* (w1Prev(eVec+1) - 2 * w1Prev(eVec) + w1Prev(eVec-1)));
    energy1w(n) = kinEnergy1w(n) + potEnergy1w(n);
    
    % Energy Mass
    kinEnergy2w(n) = M / 2 * (1/k * (w2 - w2Prev)).^2;
    potEnergy2w(n) = M / 2 * omega0^2 * ((w2 - offset) * (w2Prev - offset));
    energy2w(n) = kinEnergy2w(n) + potEnergy2w(n);
    
    % Energy Plate
    kinEnergy3w(n) = ((rhoP * H) / 2) * hP^2 * sum(sum(1/k^2 * (w3 - w3Prev).^2));
    potEnergy3w(n) = D / (2 * hP^2) * sum((Denergy * w3) .* (Denergy * w3Prev));
    energy3w(n) = kinEnergy3w(n) + potEnergy3w(n);
    
    % Connection and Collision energies
    colEnergy1w(n) = K1 / (alpha1 + 1) * 0.5 * (abs(etaw1)^(alpha1+1) + abs(etawPrev1)^(alpha1+1));
    colEnergy2w(n) = K2 / (alpha2 + 1) * 0.5 * (subplus(etaw2)^(alpha2+1) + subplus(etawPrev2)^(alpha2+1));
    
    % Total Energy
    totEnergyw(n) = energy1w(n) + energy2w(n) + energy3w(n) + colEnergy1w(n) + colEnergy2w(n);
    
    %% Update states
    psi1PrevU = psi1U; 
    psi2PrevU = psi2U; 
    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    u3Prev = u3;
    u3 = u3Next;
    
    w1Prev = w1;
    w1 = w1Next;
    
    w2Prev = w2;
    w2 = w2Next;
    
    w3Prev = w3;
    w3 = w3Next;

    quPrev = qu;
    
    out(n) = u1Next(outputPos);
    out2(n) = u3Next(outputPosPlate);
    out3(n) = u2Next;
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
       % Draw States of... 
        subplot(4,2,1);
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
%             legend(["String", "Bridge", "Body"])

        % Draw States of... 
        subplot(4,2,2);
        cla
        hold on;
        %...string
        plot(w1Next, 'Linewidth', 2);    
        %...mass
        scatter(cL, w2Next, 400, '.');    
        %...barrier
        plot([cL-1, cL+1], [w3Next(brP) w3Next(brP)], 'Linewidth', 5);  

        % Extra functions
%             ylim([-amp / 2, amp / 2]); % Set y-limit to the amplitude of the raised cosine
        grid on; 
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        title("State of the system")
%             legend(["String", "Bridge", "Body"])

        subplot(4,1,2)
        plot(u1Next-w1Next)
        %...plate
        subplot(4,2,5)
        imagesc(reshape(u3Next, [Ny,Nx])')
        xticks([])
        yticks([])
        title("Body (plate)")
        set(gca,'Fontsize', 16)

        %...plate
        subplot(4,2,6)
        imagesc(reshape(w3Next, [Ny,Nx])')
        xticks([])
        yticks([])
        title("Body (plate)")
        set(gca,'Fontsize', 16)

        subplot(4,2,7);
%         hold off
%         plot(energy1u(1:n)-energy1u(1))
%         hold on
%         plot(energy2u(1:n)-energy2u(1))
%         plot(energy3u(1:n)-energy3u(1))
%         plot(colEnergy1u(1:n) - colEnergy1u(1))
%         plot(colEnergy2u(1:n))
        % Draw Normalised energy
        plot(totEnergyu(1:n) / totEnergyu(1) - 1);
%         plot(energy3u(10:n))
%         plot(energy2u(1:n))
        title("Normalised Energy")
        
        subplot(4,2,8);

        % Draw Normalised energy
%         hold off
%         plot(energy1w(1:n)-energy1w(1))
%         hold on
%         plot(energy2w(1:n)-energy2w(1))
%         plot(energy3w(1:n)-energy3w(1))
% %         hold on
%         plot(colEnergy1w(1:n) - colEnergy1w(1))
%         plot(colEnergy2w(1:n))
        plot(totEnergyw(1:n)/totEnergyw(1) - 1);
%         plot(energy3w(10:n));
%         legend("1", "2", "3", "col1", "col2")
        title("Normalised Energy")

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