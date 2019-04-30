    clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = false;
drawSpeed = 1;
lengthSound = 1 * fs;
drawStart = 0;

%% String Variables
f0 = 110.00;    
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rho * A;  % Tension
k = 1 / fs;         % Time step
E = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappa = sqrt (E*I / (rho*A));   % Stiffness coefficient

% Damping coefficients
s0 = 2 * rho * A;
s1 = 0;
scaleFac = rho * A; % (scaling with mass/unit length)

% Grid spacing
h = sqrt((c^2*k^2 + 4*s1*scaleFac*k + sqrt((c^2*k^2 + 4*s1*scaleFac*k)^2+16*kappa^2*k^2)) / 2);
N = floor(1/h);
h = 1/N;
offset = 0.1;
u1 = zeros(N, 1) + offset;
u1Prev = zeros(N, 1) + offset;
u1Next = zeros(N, 1) + offset;

%% Mass Variables
f1 = 500;            % fundamental frequency [Hz]
w1 = 2 * pi * f1;   % angular frequency
M = 0.1;
u2 = -0.1;
u2Prev = -0.1;
% u2Next = -0.1;

%% Collision Variables
alpha = 1.3;
K = 0 * 10^6;
cL = floor (N / 6); % collision location

%% Excitation (raised cosine)
width = 20;
loc = 1/2;
startIdx = floor(floor(loc * N) - width / 2);
endIdx = floor(floor(loc * N) + width / 2);
u1(startIdx : endIdx) = u1(startIdx : endIdx) + (1 - cos(2 * pi * [0:width]' / width)) / 2;
u1Prev = u1;

%% Initialise
etaPrev = u2 - u1(cL);
psiPrev = 0;
eta = u2Prev - u1Prev(cL);
etaNext = u2 - u1(cL);

%% Initialise Energy Vectors
kinEnergy1 = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);
kinEnergy2 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);
colEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);
rOCkinEnergy1 = zeros(lengthSound, 1);
rOCpotEnergy1 = zeros(lengthSound, 1);
rOCcolEnergy1 = zeros(lengthSound, 1);
rOCdamp0Energy = zeros(lengthSound, 1);
rOCdamp1Energy = zeros(lengthSound, 1);
rOCTotEnergy1 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCcolEnergy2 = zeros(lengthSound, 1);
rOCTotEnergy2 = zeros(lengthSound, 1);

rOCTotEnergy = zeros(lengthSound, 1);

Adiv = zeros(N, 1);
v = zeros(N, 1);

vec = 3:N-2;
eVec = 2:N-1; % energy vector
J = zeros(N,1);
J(cL) = 1 / h;
outputPos = floor(N / 3);
scaleH = 1;
for n = 2:lengthSound
    %% Calculate energy of the system
    kinEnergy1(n) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h * sum((u1(3:N) - u1(2:N-1)) .* (u1Prev(3:N) - u1Prev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    kinEnergy2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2(n) = M / 2 * w1^2 * u2 * u2Prev;
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    colEnergy(n) = 1/2 * psiPrev^2;
    totEnergy(n) = energy1(n) + energy2(n) + colEnergy(n);
    
    etaPrev = eta;
    eta = etaNext; 
    
    %% Calculate g
    if alpha == 1
        g = 0;
        if eta > 0
            g = ones(N,1) * sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(eta)^((alpha - 1)/2);
    end
    %% Matrix Form
    Amat = [rho*A/k^2+g.^2/(4*h)+s0/k, -g.^2/(4*h);...
            -g.^2/4, M/k^2+g.^2/4];
    answ = [rho*A/k^2 * (2 * u1(cL) - u1Prev(cL)) + T / h^2 * (u1(cL+1) - 2 * u1(cL) + u1(cL-1)) ...
            - E * I / h^4 * (u1(cL+2) - 4 * u1(cL+1) + 6 * u1(cL) - 4 * u1(cL-1) + u1(cL-2))...
            + s0 / k * u1Prev(cL)...
            - (g.^2 / 4 * etaPrev - psiPrev .* g) / h;...
            M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2 + g.^2/4 * etaPrev - psiPrev .* g];
    solut = Amat\answ;
    etaNext = solut(2) - solut(1);
    u2Next = solut(2);
    %% Non-Matrix Form (WITHOUT DAMPING)
%     etaNextDiv =  1 + g^2 / 4 * k^2 / M + g^2 / 4 * k^2 / (rho * A * h);
%     etaNextV = (M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2 + g^2 / 4 * etaPrev - psiPrev * g) * k^2 / M...
%         - (rho * A / k^2 * (2 * u1(cL) - u1Prev(cL)) + T / h^2 * (u1(cL+1) - 2 * u1(cL) + u1(cL-1))...
%         - E * I / h^4 * (u1(cL+2) - 4 * u1(cL+1) + 6 * u1(cL) - 4 * u1(cL-1) + u1(cL-2))...
%         - (g^2 / 4 * etaPrev - psiPrev * g) / h) * k^2 / (rho * A);
%     etaNext = etaNextV / etaNextDiv;
    
    %% Update FDS
%     Adiv(vec) = rho * A / k^2;
    u1Next(vec) = (rho * A / k^2 * (2 * u1(vec) - u1Prev(vec)) + T / h^2 * (u1(vec+1) - 2 * u1(vec) + u1(vec-1)) ...
         - E * I / h^4 * (u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2))...
         + s0 / k * u1Prev(vec)...
         + J(vec) * (g^2/4 * (etaNext - etaPrev) + psiPrev * g)) / (rho * A / k^2 + s0/k);

%     u2Next = (M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2 - g^2 / 4 * (etaNext - etaPrev) - psiPrev * g) * k^2 / M;
    if drawThings
        etaNext - (u2Next - u1Next(cL))
%         u2Next - solut(2)
    end
    %% Update Psi
    psi = psiPrev + 0.5 * g .* (etaNext - etaPrev);
    
%     %% Calculate rate-of-changes in energy for checking damping (inner product of delta t-dot with scheme)     
    rOCkinEnergy1(n) = h * rho * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = h * T / (2*k*h^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
         - h * E * I / (2 * k * h^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0Energy(n) = -2 * s0 * h / (4 * k^2) * sum((u1Next - u1Prev).^2);
%     rOCdamp1Energy(n) = scaleFac * 2 * h * s1 / (2 * k^2 * h^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    energy1(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0Energy(n);
    rOCcolEnergy1(n) = 1 / (4 * k) * g * (psi + psiPrev) * (u1Next(cL) - u1Prev(cL));
    
    rOCkinEnergy2(n) = M / (2*k^3) * (u2Next - u2Prev) * (u2Next - 2 * u2 + u2Prev);
    rOCpotEnergy2(n) = -M * w1^2 / (2*k) * (u2Next - u2Prev) * u2;
    energy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n);
    
    rOCcolEnergy2(n) = -1 / (4 * k) * g * (psi + psiPrev) * (u2Next - u2Prev);
    rOCTotEnergy(n) = energy1(n) + energy2(n) - rOCcolEnergy1(n) - rOCcolEnergy2(n); %including damping so should be 0
    %% Update states
    psiPrev = psi; 
    u1Prev = u1;
    u1 = u1Next;
    u2Prev = u2;
    u2 = u2Next;
    
    out(n) = u1Next(outputPos);
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(2,1,1);
        cla
        plot(u1Next);
        hold on;
        scatter(cL, u2Next);
       
        subplot(2,1,2);
        cla
        if s0 == 0 && s1 == 0
             plot(totEnergy(10:n) / totEnergy(10) - 1);
            title("Normalised total energy");
        else
          plot(rOCTotEnergy(10:n));
%             plot(r OCTotEnergy(10:n));
%             hold on;
%             plot(rOCcolEnergy1(10:n));
%             plot(rOCcolEnergy2(10:n));
%             title("Rate of change total energy");
        end
        drawnow
    end
end
if ~drawThings
    plot(out)
end