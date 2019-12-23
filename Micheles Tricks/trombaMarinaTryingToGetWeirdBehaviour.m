%{
 The "weird behaviour" this script attempts to show is the fact that the
string gets "stuck" behind the mass in the explicit case. In terms of the
algorithm, this means that psi is negative for a prolonged period of time
In this case, this happens from around sample 2800 with a sample rate of
44100*10. 
%}

clear all;
close all;

%% Sample rate
fsScalar = 10; % 1: no weird behaviour, 10: weird behaviour, 20, no weird behaviour
fs = 44100 * fsScalar; 
k = 1/fs;

%% Excitation (cos or bowed)
exc = "cos";
plotWeirdBehaviour = true;
stopBowSample = 120000 * fs / 44100;

%% Drawing Functions
drawThings = true;
drawSpeed = 1 * fs / 44100;
lengthSound = fs * 2;
drawStart = 0;

%% String Variables
f0 = 100;
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
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

u1Next = zeros(NS, 1);
w1Next = zeros(NS, 1);
u1 = zeros(NS, 1);
w1 = zeros(NS, 1);
u1Prev = zeros(NS, 1);
w1Prev = zeros(NS, 1);

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
tol = 1e-4; % bowing NR tolerance 

%% Mass Variables
f1 = 1000;    % fundamental frequency [Hz] (< 1 / (k * pi) (< 14,037 Hz))
omega0 = 2 * pi * f1;   % angular frequency (< 2 / k (< 88,200 rad/s))
M = 0.001;

% Initialise the mass state with a velocity to make sure that 
% etaNext - etaPrev ~= 0 (which is what we are dividing by in the NR algorithm)

u2 = -0.99999e-6;
u2Prev = -1e-6;
u2Next = 0;

w2 = -0.99999e-6;
w2Prev = -1e-6;
w2Next = 0;

% Mass location along the string
massLoc = 0.8;

%% Collision Variables
cL = floor (NS * massLoc) + 1; % bridge location
alpha1 = 1.5;
K = 5 * 10^10;

%% Excitation
amp = -1e-4;

if strcmp(exc, "cos")
    width = floor(NS / 5);
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
colTol = 1e-15; % collision NR tolerance

etauPrev = u2Prev - u1Prev(cL);
etaNextPrev = etauPrev;

etaw = w2 - w1(cL);
etawPrev = w2Prev - w1Prev(cL);
etawNext = etaw;

strVec = 3:NS-2; 
for n = 1:lengthSound 
    Fb = FbInit;
    
    %% Explicit method %%
    
    % calculate eta^n and eta^{n-1}
    etau = u2 - u1(cL);
    etauPrev = u2Prev - u1Prev(cL);
     
    %% Update FDSs without connection-force and collision terms 
    u1Next(strVec) = BS(strVec, :) * u1 + CS(strVec, :) * u1Prev;
    u2Next = 2 * u2 - u2Prev - k^2 * omega0^2 * u2;    
    
    %% Calculate g^n
    g = 0;
    if alpha1 == 1
        if etau > 0
            g = sqrt(K * (alpha1+1) / 2);
        end
    else
        g = sqrt(K * (alpha1+1) / 2) * subplus(etau)^((alpha1 - 1)/2);
    end
    gSave(n) = g;

    %% Bowing
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
    
    % Stop bowing after "stopBowSample" samples
    if n < stopBowSample
        u1Next = u1Next - k^2 / (rho * A) * (Jbow * Fb * BM * qu * exp(-a*qu^2));
    end
    
    %% Solve system
    Adiv = 1 + g^2*k^2 / (4*M) + g^2*k^2 / (4*rho*A*hS);
    v = (u2Next - u2Prev - u1Next(cL) + u1Prev(cL)) / (2*k) - k*psiPrev*g / (2*M) - psiPrev*g*k / (2*rho*A*hS);
    solut = v / Adiv;
    
    psi = g*k*solut + psiPrev;
    
    %% Calculate force and add to the string and mass-spring
    forceEXP(n)  = 0.5 * (psi + psiPrev) * g;
    u1Next(cL) = u1Next(cL) + k^2 / (rho * A * hS) * forceEXP(n);
    u2Next = u2Next - k^2 / M * forceEXP(n);
    
    
    
    
    %% Newton-Raphson %%
    
    %% Update FDSs without connection-force and collision terms 
    w1Next(strVec) = BS(strVec, :) * w1 + CS(strVec, :) * w1Prev;
    w2Next = 2 * w2 - w2Prev - k^2 * omega0^2 * w2;

    %% Bowing
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
    
    % Stop bowing after "stopBowSample" samples
    if n < stopBowSample
        w1Next = w1Next - k^2 / (rho * A) * (Jbow * Fb * BM * qw * exp(-a*qw^2));
    end
       
    % caluclate eta^n and eta^{n-1}
    etaw = w2 - w1(cL);
    etawPrev = w2Prev - w1Prev(cL);
    
    %% NR collision
    b = -2 * etaw + etawPrev + k^2 * omega0^2 * w2 + T * k^2 / (rho * A * hS^2) * (w1(cL+1) - 2 * w1(cL) + w1(cL-1))...
        - E * I * k^2 / (rho * A * hS^4) * (w1(cL+2) - 4 * w1(cL+1) + 6 * w1(cL) - 4 * w1(cL-1) + w1(cL-2));
    eps = 1;
    ii = 0;

    while eps > colTol && ii < 100
        % terms to use for the quotient rule
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

    %% Calculate force and add to the string and mass-spring
    forceNR(n) = K / (alpha1 + 1) * (subplus(etawNext)^(alpha1+1) - subplus(etawPrev)^(alpha1+1)) / (etawNext - etawPrev);
    w1Next(cL) = w1Next(cL) + k^2 / (rho * A * hS) * forceNR(n);
    w2Next = w2Next - k^2 / M * forceNR(n);

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
        + E * I / 2 * 1/hS^3 * sum((w1(eVec+1) - 2 * w1(eVec) + w1(eVec-1)) ...
        .* (w1Prev(eVec+1) - 2 * w1Prev(eVec) + w1Prev(eVec-1)));
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
    
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        
        %% Draw the state of the explicitly calc
        subplot(3,3,1)
        cla
        hold on;
        %...string
        plot(u1Next, 'Linewidth', 2);    
        %...mass
        scatter(cL, u2Next, 400, '.');    
        
        grid on; 
        set(gca, 'Linewidth', 2)
        set(gca, 'Fontsize', 16)
        title("System state (explicit)")
        
        subplot(3,3,2)
        cla
        hold on;
        %...string
        plot(w1Next, 'Linewidth', 2);    
        %...mass
        scatter(cL, w2Next, 400, '.');    
        
        grid on; 
        set(gca, 'Linewidth', 2)
        set(gca, 'Fontsize', 16)
        title("System state (NR)") 
        
        subplot(3,3,3)
        cla
        hold on;
        %...string
        plot(u1Next, 'Linewidth', 2);        
        plot(w1Next, 'Linewidth', 2);    
        
        %...mass
        scatter(cL, u2Next, 400, '.');    
        scatter(cL, w2Next, 400, '.');
        
        grid on; 
        set(gca, 'Linewidth', 2)
        set(gca, 'Fontsize', 16)
        title("Both states")
        
        
        subplot(3,3,4)
        plot(forceEXP(1:n));
        title("$(\mu_{t+}\psi^{n-1/2})g^n$", 'interpreter', 'latex')
        set(gca, 'Fontsize', 16)

        subplot(3,3,5)
        plot(forceNR(1:n));
        set(gca, 'Fontsize', 16)
        title("$\delta_{t\cdot}\phi^{n} / \delta_{t\cdot}\eta^n$", 'interpreter', 'latex') 

        subplot(3,3,6)
        cla
        plot(forceEXP(2:n-1));
        hold on
        plot(forceNR (2:n-1));
        set(gca, 'Fontsize', 16)
        title("Both") 

        
        subplot(3,3,7)
        plot(totEnergyu(1:n) / totEnergyu(1) - 1);
        title("Energy u")
        set(gca, 'Fontsize', 16)

        subplot(3,3,8)
        plot(totEnergyw(1:n) / totEnergyw(1) - 1);
        set(gca, 'Fontsize', 16)
        title("Energy w") 

        subplot(3,3,9)
        cla
        plot(totEnergyu(1:n));
        hold on
        plot(totEnergyw(1:n));
        set(gca, 'Fontsize', 16)
        title("Both") 

        drawnow;
    end
    
end