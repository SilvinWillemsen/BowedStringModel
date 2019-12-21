clear all;
close all;

%% Sample rate
fs = 44100*10;
k = 1/fs;

%% Excitation (cos or bowed)
exc = "bowed";
contCosRate = 100;
plotWeirdBehaviour = true;
stopBowSample = 120000 * fs / 44100;

%% Drawing Functions
drawThings = true;
drawSpeed = 10 * fs / 44100;
lengthSound = fs * 2;
drawStart = 0;%250 * fs / 44100;
damping = true;
dampTest = false;
onlyString = false;

%% Bridge offset and location
bridgeLoc = 0.8;

%% String Variables
f0 = 100;
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
% T = 10;  % Tension
T = c^2 * rho * A;
k = 1 / fs;         % Time step
E = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappaS = sqrt (E*I / (rho*A));   % Stiffness coefficient

% Damping coefficients
s0 = 0.0;
s1 = 0;

[BS, CS, NS, hS, Dxx, Dxxxx, s0, s1, bB, bC] = unscaledCreateStringNR (rho, A, T, E, I, L, s0, s1, k);

u1 = zeros(NS, 1);
u1Prev = zeros(NS, 1);
u1Next = zeros(NS, 1);

w1 = zeros(NS, 1);
w1Prev = zeros(NS, 1);
w1Next = zeros(NS, 1);

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
colTol = 1e-10;

%% Mass Variables
f1 = 1000;    % fundamental frequency [Hz] (< 1 / (k * pi) (< 14,037 Hz))
omega0 = 2 * pi * f1;   % angular frequency (< 2 / k (< 88,200 rad/s))
M = 0.001;

u2 = -0.99999e-6;
u2Prev = -1e-6;
u2Next = 0;

w2 = -0.99999e-6;
w2Prev = -1e-6;
w2Next = 0;

%% Collision Variables
cL = floor (NS * bridgeLoc) + 1; % bridge location
alpha1 = 1.5;
K = 5 * 10^10;

%% Excitation
amp = -1e-4;

if strcmp(exc, "cos")
    width = 10;
    loc = 1/2;
    startIdx = floor(floor(loc * NS) - width / 2) + 1;
    endIdx = floor(floor(loc * NS) + width / 2) + 1;
    u1Next(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
    w1Next(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;

end
u1 = u1Next;
u1Prev = u1Next;
w1 = w1Next;
w1Prev = w1Next;

%% Initialise
psiPrev = 0;

%% Initialise Energy Vectors
kinEnergyu1 = zeros(lengthSound, 1);
potEnergyu1 = zeros(lengthSound, 1);
energyu1 = zeros(lengthSound, 1);

kinEnergyu2 = zeros(lengthSound, 1); 
potEnergyu2 = zeros(lengthSound, 1);
energyu2 = zeros(lengthSound, 1);

kinEnergy3 = zeros(lengthSound, 1);
potEnergy3 = zeros(lengthSound, 1);
energy3 = zeros(lengthSound, 1);

colEnergyu1 = zeros(lengthSound, 1);
colEnergy2 = zeros(lengthSound, 1);
totEnergyu = zeros(lengthSound, 1);

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
Ibr = zeros(NS,1);
Ibr(cL) = 1;
Jbr = Ibr / hS;

Ibow = zeros(NS,1);
Ibow(bP) = 1;
Jbow = Ibow / hS;

outputPos = floor(NS / 6);
prog = 0;

rampVal = 10000;
etauPrev = u2Prev - u1Prev(cL);
etaNextPrev = etauPrev;

etaw = w2 - w1(cL);
etawPrev = w2Prev - w1Prev(cL);
etawNext = etaw;

for n = 1:lengthSound 

    Fb = FbInit;
    etau = u2 - u1(cL);
    etauSave(n) = etau;
    
    %% Update FDSs without connection-force and collision terms 
    strVec = 3:NS-2;
    
    u1Next(strVec) = BS(strVec, :) * u1 + CS(strVec, :) * u1Prev;
    u2Next = 2 * u2 - u2Prev - k^2 * omega0^2 * u2;    
    
    %% EXPLICIT
    
    g = 0;
    if alpha1 == 1
        if etau > 0
            g = sqrt(K * (alpha1+1) / 2);
        end
    else
        g = sqrt(K * (alpha1+1) / 2) * subplus(etau)^((alpha1 - 1)/2);
    end
    gSave(n) = g;
    
    %% Bowing u
    if exc == "bowed"
        b = 2/k * Vb + Ibow' * bB * u1 + Ibow' * bC * u1Prev;
        eps = 1;
        ii = 0;

        while eps>tol && ii < 100
            qu=quPrev-(2*quPrev/k + Fb * BM /(rho * A * hS)*quPrev*exp(-a*quPrev^2)+b)/...
             (Fb*BM/(rho * A * hS)*(1-2*a*quPrev^2)*exp(-a*quPrev^2)+2/k);
            eps = abs(qu-quPrev);
            quPrev = qu;
            ii = ii + 1;
        end
    else
        qu = 0;
    end
%     iiSave(n) = ii;
    if n < stopBowSample
        u1Next = u1Next - k^2 / (rho * A) * (Jbow * Fb * BM * qu * exp(-a*qu^2));
    end
    Adiv = 1 + g^2*k^2 / (4*M) + g^2*k^2 / (4*rho*A*hS);
    v = (u2Next - u2Prev - u1Next(cL) + u1Prev(cL)) / (2*k) - k*psiPrev*g / (2*M) - psiPrev*g*k / (2*rho*A*hS);
     
    solut = v / Adiv;
    
    psi = g*k*solut + psiPrev;
    if g == 0
        psiSave(n) = 0;
    else
        psiSave(n) = psi;
    end
    
    forceEXP(n)  = 0.5 * (psi + psiPrev) * g;
    u1Next(cL) = u1Next(cL) + k^2 / (rho * A * hS) * forceEXP(n);
    u2Next = u2Next - k^2 / M * forceEXP(n);

    etauNext = u2Next - u1Next(cL);
    
    if etauNext - etauPrev == 0
        diffSave(n) = 0;
    else
        diffSave(n) = 1/k * (psi - psiPrev) / (1/(2*k) * (etauNext - etauPrev));
    end
    
    %% NR
    w1Next(strVec) = BS(strVec, :) * w1 + CS(strVec, :) * w1Prev;
    w2Next = 2 * w2 - w2Prev - k^2 * omega0^2 * w2;
    % Bowing w
    if exc == "bowed"
        b = 2/k * Vb + Ibow' * bB * u1 + Ibow' * bC * u1Prev;
        eps = 1;
        ii = 0;

        while eps>tol && ii < 100
            qw=qwPrev-(2*qwPrev/k + Fb * BM /(rho * A * hS)*qwPrev*exp(-a*qwPrev^2)+b)/...
             (Fb*BM/(rho * A * hS)*(1-2*a*qwPrev^2)*exp(-a*qwPrev^2)+2/k);
            eps = abs(qw-qwPrev);
            qwPrev = qw;
            ii = ii + 1;
        end
    else
        qw = 0;
    end
%     iiSave(n) = ii;
    if n < stopBowSample
        w1Next = w1Next - k^2 / (rho * A) * (Jbow * Fb * BM * qw * exp(-a*qw^2));
    end
   
    %% NR collision
    etaw = w2 - w1(cL);
    etawPrev = w2Prev - w1Prev(cL);
    
    b = -2 * etaw + etawPrev + k^2 * omega0^2 * w2 + T * k^2 / (rho * A * hS^2) * (w1(cL+1) - 2 * w1(cL) + w1(cL-1))...
        - E * I * k^2 / (rho * A * hS^4) * (w1(cL+2) - 4 * w1(cL+1) + 6 * w1(cL) - 4 * w1(cL-1) + w1(cL-2));
    eps = 1;
    ii = 0;

    while eps > colTol && ii < 100
        f = subplus(etawNext)^(alpha1+1) - subplus(etawPrev)^(alpha1+1);
        df = (alpha1+1) * subplus(etawNext)^alpha1 * 0.5 * (1+sign(etawNext));
        gg = etawNext - etawPrev;
        dgg = 1;
        
        coeff = (k^2 / (rho * A * hS) + k^2 / M) * K / (alpha1+1);
        G = etawNext + coeff * (subplus(etawNext)^(alpha1+1)-subplus(etawPrev)^(alpha1+1)) / (etawNext - etawPrev) + b;
        Gdiff = 1 + coeff * (df * gg - f * dgg) / gg^2;
        etaNextNew = etawNext - G / Gdiff;
        
        eps = abs(etaNextNew - etawNext);
        etawNext = etaNextNew;
        ii = ii + 1;
    end

    
    forceNR(n) = K / (alpha1 + 1) * (subplus(etawNext)^(alpha1+1) - subplus(etawPrev)^(alpha1+1)) / (etawNext - etawPrev);
%     if collEffect ~= 0
%         disp("wait");
%     end
    w1Next(cL) = w1Next(cL) + k^2 / (rho * A * hS) * forceNR(n);
    w2Next = w2Next - k^2 / M * forceNR(n);
%     disp(ii + " " + (etawNext - (w2Next - w1Next(cL))))

    if drawThings && drawStart == 0
%         disp("Difference between etaNext and u3Next(brP) - u2Next (should be 0)")
%         etaNext - (u3Next(brP) - u2Next)
    elseif mod(n,floor(lengthSound / 100)) == 0
        prog = prog + 1;
        disp("Progress: " + prog + "%")
    end
    
    %% Calculate energy of u
    % Energy String
    kinEnergyu1(n) = rho * A / 2 * hS * sum((1/k * (u1 - u1Prev)).^2);
    potEnergyu1(n) = T / 2 * 1/hS * sum((u1(3:NS) - u1(2:NS-1)) .* (u1Prev(3:NS) - u1Prev(2:NS-1)))...
        + E * I / 2 * 1/hS^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energyu1(n) = kinEnergyu1(n) + potEnergyu1(n);
    
    % Energy Mass
    kinEnergyu2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergyu2(n) = M / 2 * omega0^2 * (u2 * u2Prev);
    energyu2(n) = kinEnergyu2(n) + potEnergyu2(n);
    
    % Connection and Collision energies
    colEnergyu1(n) = psiPrev^2 / 2;
    
    % Total Energy
    totEnergyu(n) = energyu1(n) + energyu2(n) + colEnergyu1(n);
    
    %% Calculate energy of w
    % Energy String
    kinEnergyw1(n) = rho * A / 2 * hS * sum((1/k * (w1 - w1Prev)).^2);
    potEnergyw1(n) = T / 2 * 1/hS * sum((w1(3:NS) - w1(2:NS-1)) .* (w1Prev(3:NS) - w1Prev(2:NS-1)))...
        + E * I / 2 * 1/hS^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energyw1(n) = kinEnergyw1(n) + potEnergyw1(n);
    
    % Energy Mass
    kinEnergyw2(n) = M / 2 * (1/k * (w2 - w2Prev)).^2;
    potEnergyw2(n) = M / 2 * omega0^2 * (w2 * w2Prev);
    energyw2(n) = kinEnergyw2(n) + potEnergyw2(n);
    
    % Connection and Collision energies
    colEnergyw1(n) = K / (alpha1 + 1) * 0.5 * (subplus(etaw)^(alpha1+1) + subplus(etawPrev)^(alpha1+1));
    
    % Total Energy
    totEnergyw(n) = energyw1(n) + energyw2(n) + colEnergyw1(n);
    
    %% Update states
    psiPrev = psi; 
    etauPrev = etau;

    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    w1Prev = w1;
    w1 = w1Next;
    
    w2Prev = w2;
    w2 = w2Next;

    quPrev = qu;
    qwPrev = qw;
    
    out(n) = u1Next(outputPos);
    out3(n) = u2Next;
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(3,1,1)
        cla
        hold on;
        %...string
        plot(u1Next, 'Linewidth', 2);    
        plot(w1Next, 'Linewidth', 2);    
        %...mass
        scatter(cL, u2Next, 400, '.');    
        scatter(cL, w2Next, 400, '.');    
        % Extra functions
%             ylim([-1e-4, 1e-4]); % Set y-limit to the amplitude of the raised cosine
        grid on; 
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        title("State of the system")
%             legend(["String", "Bridge", "Body"])
        subplot(3,1,2)
%         title("Normalised energy u")
%         plot(totEnergyu(1:n) / totEnergyu(1) - 1);
        plot(forceEXP(1:n))
        subplot(3,1,3)
%         title("Normalised energy w")
%         plot(totEnergyw(1:n) / totEnergyw(1) - 1);
        plot(forceNR(1:n));
    end
    drawnow;
end
posOut1 = out - min(out);
totOut1 = (posOut1/max(abs(posOut1)) - 0.5) * 2;
posOut2 = out2 - min(out2);
totOut2 = (posOut2/max(abs(posOut2)) - 0.5) * 2;
posOut3 = out3 - min(out3);
totOut3 = (posOut3/max(abs(posOut3)) - 0.5) * 2;
totOut4 = (gSave/max(abs(gSave)) - 0.5) * 2;
totOut = totOut1/5 + totOut2 + totOut3;

if ~drawThings
    if onlyString
        plot(out);
    else
        plot(totOut);
    end
end
plot(g2Save)