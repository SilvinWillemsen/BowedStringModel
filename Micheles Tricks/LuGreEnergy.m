clear all;
close all;

drawString = true;
drawSpeed = 10;
exc = "bowed";
fs = 44100;
k = 1 / fs;     % Time step

%% Extra Settings
lengthSound = fs * 2;

%% String variables
L = 1;
f0 = 196;    
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rho * A;  % Tension
E = 2e11;
In = r^4 * pi / 4;
s0 = 0.5;
s1 = 0.005;
[B, C, N, h, Dxx, Dxxxx, s0, s1, bB, bC] = unscaledCreateStringNR (rho, A, T, E, In, L, s0, s1, k);

kappa = sqrt(E * In / (rho * A));

courantNoS = c^2 * k^2 / h^2 + 4 * kappa^2 * k^2 / h^4

u = zeros(N, 1);

%% Raised cosine
if exc == "cos"
    amp = 0.1;
    width = 10;
    loc = 2/5;
    startIdx = floor(floor(loc * N) - width / 2);
    endIdx = floor(floor(loc * N) + width / 2);
    u(startIdx : endIdx) = u(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
end
uNext = zeros(N, 1);
uPrev = u;

u2 = zeros(N, 1);

%% Raised cosine
if exc == "cos"
    amp = 0.1;
    width = 10;
    loc = 2/5;
    startIdx = floor(floor(loc * N) - width / 2);
    endIdx = floor(floor(loc * N) + width / 2);
    u2(startIdx : endIdx) = u2(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
end 
u2Next = zeros(N, 1);
u2Prev = u2;

VbInit = 0.2;
if exc == "cos"
    Fb = 0.00;
else
    Fb = 1;
end
a = 100;

BM = sqrt(2 * a) * exp(1/2);
BP = 1/pi;      % Bowing Position
bp = floor(BP * N);

% qPrev = -Vb;
qPrev = -VbInit;
qPrev2 = -VbInit;
I = zeros(N,1);
I(bp) = 1;
rampVal = 10000;
out = zeros(lengthSound, 1);
J = 1/h * I;

%% LuGre model parameters
muS=0.8;               % static friction coeff
muC=0.3;              % dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
strv=1e-1;             % "stribeck" velocity

Fn = 1;
fS=muS*Fn;             % stiction force
fC=muC*Fn;             % coulomb force

v = 0.5;
espon = exp(-(v/strv)^2);


sig0 = 1e4;
gv = 1/sig0 * (fC + (fS - fC) * espon);

% sig1 = 4 * sig0 * gv / abs(v);
% sig1 = 0.01 * sqrt(sig0);
sig1 = 0;
sig2 = 0.0;

% Init NR 
v = 0;
z = 0;
zdot = 0;
zPrev = 0;
anPrev = 0;

tol = 1e-9;

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
potEnergyBow = zeros(lengthSound, 1);
energy = zeros(lengthSound, 1);

rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCdamp0StringEnergy = zeros(lengthSound, 1);
rOCdamp1StringEnergy = zeros(lengthSound, 1);

rOCbowEnergy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);
rOCbowEnergyInput = zeros(lengthSound, 1);
rOCbristleEnergy = zeros(lengthSound, 1);
rOCbowEnergy = zeros(lengthSound, 1);
rOCdispEnergy = zeros(lengthSound, 1);

for n = 1 : lengthSound
    Vb = VbInit;
    
    % Newton-Raphson
    b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
    eps = 1;
    i = 0;
    while eps>tol && i < 100
        espon = exp(-(v/strv)^2);
        gv = 1/sig0 * (fC + (fS - fC) * espon);
        
        zdot = v - (abs(v) / gv) * z;
        an = 2/k * (z - zPrev) - anPrev;
        
        g1 = (2/k + 2 * s0) * v + (sig0 * z + sig1 * zdot + sig2 * v) / (rho * A * h) + b;
        g2 = zdot - an;
        
        %% Compute derivatives
        dgv = -2 * v / (sig0 * strv^2) * (fS - fC) * espon;
        
        dzdot_z = -abs(v) / gv;
        dzdot_v = 1 - z * (sign(v) * gv - dgv * abs(v)) / gv^2;
        
        dg1v = 2/k + 2 * s0 + sig1 * dzdot_v / (rho * A * h) + sig2 / (rho * A * h);
        dg1z = sig0 / (rho * A * h) + sig1 * dzdot_z / (rho * A * h);
        dg2v = dzdot_v;
        dg2z = dzdot_z - 2 / k;
        
        Jac = [dg1v, dg1z; ...
            dg2v, dg2z];
        
        prevSolut = [v; z];
        solut = [v; z] - Jac \ [g1; g2];

        v = solut(1);
        z = solut(2);
        
        eps = norm(solut - prevSolut);
        
        i = i + 1;
        
    end
    vSave(n) = v;
    zSave(n) = z;
    iSave(n) = i;
    
    F = sig0 * z + sig1 * zdot + sig2 * v;

    excitation = J * F;
    uNext = B * u + C * uPrev;
    
    if exc == "bowed"
        uNext = uNext - excitation / (rho * A / k^2 + s0 / k);
    end
    
    vrelSave(n) = (uNext(bp) - uPrev(bp)) / (2*k) - Vb;
%     vSave(n) - vrelSave(n)
    eVec = 2:N-1;
    vec = 3:N-2;
    
    kinEnergy(n) = rho * A * h / 2 * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2*h) * sum((u(eVec + 1) - u(eVec)) .* (uPrev(eVec + 1) - uPrev(eVec))) + ...
        E * In / (2 * h^3) * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) ...
        .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)));
    potEnergyBow(n) = sig0 / 2 * z^2;
    energy(n) = kinEnergy(n) + potEnergy(n);
    
    rOCkinEnergy(n) = rho * A * h / (2 * k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = T / (2 * k * h) * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) .* (uNext(eVec) - uPrev(eVec)))...
         - h * E * In / (2 * k * h^4) * sum((u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2)) .* (uNext(vec) - uPrev(vec)));
    rOCdamp0StringEnergy(n) = -2 * s0 * h / (4 * k^2) * sum((uNext - uPrev).*(uNext - uPrev));
    rOCdamp1StringEnergy(n) = 2 * h * s1 / (2 * k^2 * h^2) * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1) - uPrev(eVec+1) + 2 * uPrev(eVec) - uPrev(eVec-1)) .* (uNext(eVec) - uPrev(eVec)));
    rOCbowEnergyInput(n) = F * Vb;
    rOCbristleEnergy(n) = sig0 * z * zdot;
    %     rOCbowEnergy(n) = -F * (uNext(bp) - uPrev(bp)) / (2 * k);
    rOCdispEnergy(n) = sig1 * (zdot + 1/2 * z * abs(v) / gv)^2 + z^2 * abs(v) / gv * (sig0 - sig1 * abs(v)/(4*gv)) + sig2 * v^2;
    if exc == "bowed"
        rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCdamp0StringEnergy(n) - rOCdamp1StringEnergy(n) + rOCbowEnergyInput(n) + rOCbristleEnergy(n) + rOCdispEnergy(n);
    else
        rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCdamp0StringEnergy(n) - rOCdamp1StringEnergy(n);
    end
%     rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCdamp0StringEnergy(n) - rOCdamp1StringEnergy(n) - rOCbowEnergy(n);

    if mod(n,drawSpeed) == 0 && drawString == true %&& n > 10000
        subplot(4,1,1);
        plot(uNext);
        
%         subplot(3,2,2);
%         plot(u2Next);
%         
        subplot(4,1,2);
        plot(energy(10:n)/energy(10) - 1);
        
%         subplot(3,2,4);
%         plot(energy2(10:n) / energy2(10) - 1);
%         
        subplot(4,1,3);
        cla;
        plot(rOCkinEnergy(10:n) - rOCpotEnergy(10:n) - rOCdamp0StringEnergy(10:n) - rOCdamp1StringEnergy(10:n));
        hold on;
        plot(rOCbowEnergyInput(10:n));
        plot(rOCdispEnergy(10:n));
        plot(rOCbristleEnergy(10:n));
%         plot(rOCtotEnergy(10:n))      
        subplot(4,1,4);
        plot(rOCtotEnergy(10:n))  
%         subplot(3,2,6);
%         plot(rOCtotEnergy2(10:n))
        drawnow;
    end
    
    uPrev = u;
    anPrev = an;
    zPrev = z;
    u = uNext;
    
    out(n) = u(floor(2*N / 3));
end

% totEnergy = kinEnergy + potEnergy;
% totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
% maxTotEnergy = max(totEnergy) - min(totEnergy);
% % plot (totEnergy);
plot(out);