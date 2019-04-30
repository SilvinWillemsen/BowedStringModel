clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Drawing Functions
drawThings = true;
drawSpeed = 1;
lengthSound = 2 * fs;
drawStart = 0;

%% String Variables
f0 = 200;    
rho = 7850;
r = 0.0009;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
T = c^2 * rho * A;  % Tension
k = 1 / fs;         % Time step
E = 2e11;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappa = sqrt (E*I / (rho*A));   % Stiffness coefficient
 
% Damping coefficients
s0 = 0.1;
s1 = 0.05;

% Grid spacing
h = sqrt((c^2*k^2 + 4*s1*k + sqrt((c^2*k^2 + 4*s1*k)^2+16*kappa^2*k^2)) / 2);
N = floor(1/h);
h = 1/N;
offset = 0.0;
u1 = zeros(N, 1) + offset;
u1Prev = zeros(N, 1) + offset;
u1Next = zeros(N, 1) + offset;

lambdaSq = c^2*k^2/h^2;
muSq = kappa^2 * k^2 / h^4;

courantNo = lambdaSq + 4 * muSq
%% Bridge offset
b = 0.5;
barrier = -2e-6;

%% Mass Variables
f1 = 0;            % fundamental frequency [Hz]
w1 = 2 * pi * f1;   % angular frequency
M = 0.001;
u2 = barrier / 2;
u2Prev = barrier / 2;
% u2Next = -0.1;

%% Collision Variables
alphaMS = 1.3;
KMS = 5 * 10^6;
cL = floor (N / 8); % bridge location

alphaMB = 1.3;
KMB = 5 * 10^8;

%% Initial condition string because of bridge
% u1(3:cL) = (1:cL-2) / (cL - 2) * (b+2*eps);
% u1(cL:N-1) = (N-1-cL:-1:0) / (N-1-cL) * (b+2*eps);
% u1Prev = u1;

%% Excitation
excitation = "cos";
if excitation == "bowed"
    alphaBow = 100;
    bp = floor(3 * N / 4);
    BM = sqrt(2 * alphaBow) * exp(1/2);
    FbInit = 1;
    VbInit = 0.2;
    %variables micheles model
%     FbInit = 0.005;
%     VbInit = 0.03;
    bowRamp = 5000; % number of samples that Fb and Vb will ramp up
    if bowRamp == 0
        Vrel = (u1(bp) - u1Prev(bp)) / k - VbInit;
    else
        Vrel = 0;
    end
    model = "NR";
    tol = 1e-4;
    amp = b;
elseif excitation == "cos"
    width = 10;
    loc = 1/5;
    startIdx = floor(floor(loc * N) - width / 2);
    endIdx = floor(floor(loc * N) + width / 2);
    amp = 0.00001;
    u1(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
    u1Prev = u1;
elseif excitation == "plucked"
    u1(3:cL) = (1:cL-2) / (cL - 2) * (b+2*eps);
    u1(cL:N-1) = (N-1-cL:-1:0) / (N-1-cL) * (b+2*eps);
    u1Prev = u1;
    amp = b;
end

u1NextTest = u1Next;
u1Test = u1;
u1PrevTest = u1Prev;

%% Initialise
etaPrevMS = u2 - u1(cL);
psiPrevMS = 0;
etaMS = u2Prev - u1Prev(cL);
etaNextMS = u2 - u1(cL);

etaPrevMB = barrier - u2;
psiPrevMB = 0;
etaMB = barrier - u2Prev;
etaNextMB = barrier - u2;

%% Initialise Energy Vectors
kinEnergy1 = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);
kinEnergy2 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
energy2 = zeros(lengthSound, 1);
colEnergyMS = zeros(lengthSound, 1);
colEnergyMB = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);
rOCkinEnergy1 = zeros(lengthSound, 1);
rOCpotEnergy1 = zeros(lengthSound, 1);
rOCcolEnergyMS1 = zeros(lengthSound, 1);
rOCdamp0Energy = zeros(lengthSound, 1);
rOCdamp1Energy = zeros(lengthSound, 1);
rOCTotEnergy1 = zeros(lengthSound, 1);

rOCkinEnergy2 = zeros(lengthSound, 1);
rOCpotEnergy2 = zeros(lengthSound, 1);
rOCcolEnergyMB = zeros(lengthSound, 1);
rOCTotEnergy2 = zeros(lengthSound, 1);

rOCTotEnergy = zeros(lengthSound, 1);

Adiv = zeros(N, 1);
v = zeros(N, 1);

vec = 3:N-2;
eVec = 2:N-1; % energy vector
Jbr = zeros(N,1);
Jbr(cL) = 1 / h;


Jbow = zeros(N,1);
if strcmp(excitation, "bowed")
    Jbow(bp) = 1 / h;
end
outputPos = floor(N / 6);
scaleH = 1;
figure;
grid on;
set(gca, 'Linewidth', 2)


lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;
qPrev=0;

for n = 2:lengthSound
    %% Calculate energy of the system
    kinEnergy1(n) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h * sum((u1(3:N) - u1(2:N-1)) .* (u1Prev(3:N) - u1Prev(2:N-1)))...
        + E * I / 2 * 1/h^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    kinEnergy2(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2(n) = M / 2 * w1^2 * (u2 * u2Prev);
    energy2(n) = kinEnergy2(n) + potEnergy2(n);
    colEnergyMS(n) = 1/2 * psiPrevMS^2;
    colEnergyMB(n) = 1/2 * psiPrevMB^2;
    totEnergy(n) = energy1(n) + energy2(n) + colEnergyMS(n) + colEnergyMB(n);
    
    etaPrevMS = etaMS;
    etaMS = etaNextMS; 
    
    etaPrevMB = etaMB;
    etaMB = etaNextMB; 
    
    %% Calculate g for Mass String
    if alphaMS == 1
        gMS = 0;
        if etaMS > 0
            gMS = sqrt(KMS * (alphaMS+1) / 2);
        end
    else
        gMS = sqrt(KMS * (alphaMS+1) / 2) * subplus(etaMS)^((alphaMS - 1)/2);
    end
    
    %% Calculate g for Mass Barrier
    if alphaMB == 1
        gMB = 0;
        if etaMB > 0
            gMB = sqrt(KMB * (alphaMB+1) / 2);
        end
    else
        gMB = sqrt(KMB * (alphaMB+1) / 2) * subplus(etaMB)^((alphaMB - 1)/2);
    end
    gMBsave(n) = gMB;
    
    if strcmp(excitation, "bowed")
        %% Bow Ramp
        if n < bowRamp
            Vb = VbInit * n / bowRamp;
            Fb = FbInit * n / bowRamp;
        else
            Vb = VbInit;
            Fb = FbInit;
        end
        if strcmp(model, "NR")
            b = 2 * rho * A / k * Vb - (2 * rho * A /k^2)*(u1(bp)-u1Prev(bp))...
                - T/h^2 * (u1(bp+1) - 2*u1(bp) + u1(bp-1))...
            + E * I/h^4 * (u1(bp+2) - 4*u1(bp+1) + 6*u1(bp) - 4*u1(bp-1) + u1(bp-2)) + 2 * s0 * rho * A * Vb ...
            - (2*s1 * rho * A/(k*h^2)) * ((u1(bp+1) - 2 * u1(bp) + u1(bp-1)) ...
              - (u1Prev(bp+1) - 2 * u1Prev(bp) + u1Prev(bp-1)));
            epsilon = 1;
            i = 0;
            qPrev = 0;
            while epsilon>tol
                q=qPrev-(Fb*BM*qPrev*exp(-alphaBow*qPrev^2)+(2 * rho * A / k + 2*s0 * rho * A)*qPrev+b)/...
                 (Fb*BM*(1-2*alphaBow*qPrev^2)*exp(-alphaBow*qPrev^2)+2 * rho * A / k + 2*s0 * rho * A);
                epsilon = abs(q-qPrev);
                qPrev = q;
                i = i + 1;
                if i > 10000
                    disp("whut")
                end
            end
            qSave(n) = q;
            exciteSignal = -Jbow(vec) * Fb * BM * q * exp(-alphaBow * q^2);
        else
            %% Estimate new relative velocity
            Vrel = 2 / k * (u1(bp) - u1Prev(bp)) - Vrel - 2 * Vb;
            Btot = Fb * sqrt(2 * alphaBow) * exp(-alphaBow * Vrel^2 + 1/2);
            exciteSignal = Jbow(vec) .* (Btot / (2*k) * u1Prev(vec) + Btot * Vb);
        end
        Btot = 0;
        Vb = 0;
    else
        Btot = 0;
        exciteSignal = 0;
    end
    %% Matrix Form
    Amat = [rho*A/k^2 + gMS.^2/(4*h) + (rho * A * s0)/k,...
            -gMS.^2/(4*h);...
            -gMS.^2/4,...
            M/k^2 + gMS.^2/4 + gMB.^2/4];
    answ = [rho*A/k^2 * (2 * u1(cL) - u1Prev(cL)) + T / h^2 * (u1(cL+1) - 2 * u1(cL) + u1(cL-1)) ...
            - E * I / h^4 * (u1(cL+2) - 4 * u1(cL+1) + 6 * u1(cL) - 4 * u1(cL-1) + u1(cL-2))...
            + rho * A * s0 / k * u1Prev(cL)...
            + 2 * rho * A * s1 / (k * h^2) * (u1(cL+1) - 2 * u1(cL) + u1(cL-1) - u1Prev(cL+1) + 2 * u1Prev(cL) - u1Prev(cL-1))...
            + (-gMS.^2 / 4 * etaPrevMS + psiPrevMS .* gMS) / h;...
            M / k^2 * (2 * u2 - u2Prev) - M * w1^2 * u2 + gMS.^2/4 * etaPrevMS - psiPrevMS .* gMS...
            + gMB.^2/4 * u2Prev + psiPrevMB * gMB];
    solut = Amat\answ;
    
    etaNextMS = solut(2) - solut(1);
    etaNextMB = barrier - solut(2);
    u2Next = solut(2);
    
    %% Update FDS
    u1Next(vec) = (rho * A / k^2 * (2 * u1(vec) - u1Prev(vec)) + T / h^2 * (u1(vec+1) - 2 * u1(vec) + u1(vec-1)) ...
         - E * I / h^4 * (u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2))...
         + rho * A * s0 / k * u1Prev(vec)...
         + 2 * rho * A * s1 / (k * h^2) * (u1(vec+1) - 2 * u1(vec) + u1(vec-1) - u1Prev(vec+1) + 2 * u1Prev(vec) - u1Prev(vec-1))...
         + Jbr(vec) * (gMS^2/4 * (etaNextMS - etaPrevMS) + psiPrevMS * gMS)...
         + exciteSignal) ./ (rho * A / k^2 + rho * A * s0/k + Jbow(vec) * Btot / (2*k));
%      u1NextTest(vec) = (2 * u1Test(vec) - u1PrevTest(vec) + lambdaSq * (u1Test(vec+1) - 2 * u1Test(vec) + u1Test(vec-1))...
%          - muSq *  (u1Test(vec+2) - 4 * u1Test(vec+1) + 6 * u1Test(vec) - 4 * u1Test(vec-1) + u1Test(vec-2))...
%          + s0 * k * u1PrevTest(vec)...
%          + 2 * s1 * k / h^2 * (u1Test(vec+1) - 2 * u1Test(vec) + u1Test(vec-1) - u1PrevTest(vec+1) + 2 * u1PrevTest(vec) - u1PrevTest(vec-1))) / (1 + s0 * k);
%          
    if drawThings
        etaNextMS - (u2Next - u1Next(cL));
        etaNextMB - (b - u2Next);
    end
    %% Update Psi
    psiMS = psiPrevMS + 0.5 * gMS .* (etaNextMS - etaPrevMS);
    psiMB = psiPrevMB + 0.5 * gMB .* (etaNextMB - etaPrevMB);
    
%     %% Calculate rate-of-changes in energy for checking damping (inner product of delta t-dot with scheme)     
    rOCkinEnergy1(n) = h * rho * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy1(n) = h * T / (2*k*h^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
         - h * E * I / (2 * k * h^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0Energy(n) = -2 * rho * A * s0 * h / (4 * k^2) * sum((u1Next - u1Prev).^2);
    rOCdamp1Energy(n) = 2 * h * rho * A * s1 / (2 * k^2 * h^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1) - u1Prev(eVec+1) + 2 * u1Prev(eVec) - u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    rOCenergy1(n) = rOCkinEnergy1(n) - rOCpotEnergy1(n) - rOCdamp0Energy(n) - rOCdamp1Energy(n);
    rOCcolEnergyMS1(n) = gMS / (4 * k) * (psiMS + psiPrevMS) * (u1Next(cL) - u1Prev(cL));
    rOCcolEnergyMS2(n) = -gMS / (4 * k) * (psiMS + psiPrevMS) * (u2Next - u2Prev);
    
    rOCkinEnergy2(n) = M / (2 * k^3) * (u2Next - 2 * u2 + u2Prev) * (u2Next - u2Prev);
    rOCpotEnergy2(n) = -M * w1^2 / (2*k) * (u2Next - u2Prev) * u2;
    rOCenergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n);
    
    rOCcolEnergyMB(n) = gMB / (4 * k) * (psiMB + psiPrevMB) * (u2Next - u2Prev);
    rOCTotEnergy(n) = rOCenergy1(n) + rOCenergy2(n) - rOCcolEnergyMS1(n) - rOCcolEnergyMS2(n) - rOCcolEnergyMB(n); %including damping so should be 0
    %% Update states
    psiPrevMS = psiMS; 
    psiPrevMB = psiMB; 
    u1Prev = u1;
    u1 = u1Next;
    u1PrevTest = u1Test;
    u1Test = u1NextTest;
    u2Prev = u2;
    u2 = u2Next;
    
    out(n) = u1Next(outputPos);
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(2,1,1);
        cla
        plot(u1Next, 'Linewidth', 2);
        hold on;
%         plot(u1NextTest, 'Linewidth', 2);
        scatter(cL, u2Next, 400, '.');
        plot([cL-1, cL+1], [barrier barrier], 'Linewidth', 5);
        
        grid on; 
        set(gca, 'Linewidth', 2, 'Fontsize', 16)
        if ~strcmp(excitation, "bowed")
            ylim([-amp, amp])
        end
        legend(["String", "Bridge", "Body"])
        subplot(2,1,2);
%         cla
%         if s0 == 0 && s1 == 0
% %             plot(kinEnergy2(10:n))
% %             hold on;
%             plot(totEnergy(10:n) / totEnergy(10) - 1)
% %             plot(potEnergy2(10:n))
% %             plot(colEnergyMB(10:n))
% %             plot(totEnergy(10:n) / totEnergy(10) - 1);
% % plot(gMBsave)
% %             title("Normalised total energy");
%         else
% %             plot(rOCkinEnergy2(10:n));
% %             hold on; 
% %             plot(rOCpotEnergy2(10:n));
% %             plot(rOCdamp1Energy(10:n));
% %             hold on;
            plot(rOCTotEnergy(10:n));
% %             hold on;
% %             plot(rOCcolEnergy1(10:n));
% %             title("Rate of change total energy");
%         end
        drawnow
    end
end
if ~drawThings
    plot(out)
end