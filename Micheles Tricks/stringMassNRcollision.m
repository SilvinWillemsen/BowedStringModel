clear all;
close all;

%% Sample rate
fs = 44100;
k = 1/fs;

%% Excitation (cos or bowed)
exc = "none";
contCosRate = 100;
plotWeirdBehaviour = true;
stopBowSample = 120000 * fs / 44100;

%% collisiontype (NR or nonIt)
colType = "NR" 

%% Drawing Functions
drawThings = true;
drawSpeed = 1 * fs / 44100;
lengthSound = fs * 2;
drawStart = 0;
damping = true;
dampTest = false;
onlyString = false;

%% Bridge offset and location
bridgeLoc = 0.8;

%% String Variables
f0 = 50;
rho = 7850;
r = 0.0005;
A = r^2 * pi;
c = f0 * 2;         % Wave-Speed
% T = 10;  % Tension
T = c^2 * rho * A;
k = 1 / fs;         % Time step
E = 0;           % Young's modulus
I = r^4 * pi / 4;   % Area moment of inertia
L = 1;              % String Length
kappaS = sqrt (E*I / (rho*A));   % Stiffness coefficient

% Damping coefficients
s0 = 0.0;
s1 = 0;

[BS, CS, NS, h, Dxx, Dxxxx, s0, s1, bB, bC] = unscaledCreateStringNR (rho, A, T, E, I, L, s0, s1, k);

u1 = zeros(NS, 1);
u1Prev = zeros(NS, 1);
u1Next = zeros(NS, 1);

% e = ones(NS-1,1);
% Dxx = (1/hS^2)*spdiags([e -2*e e],-1:1,NS-1,NS-1);
% 
courantNoS = c^2 * k^2 / h^2 + 4 * kappaS^2 * k^2 / h^4

%% Bowing terms
bP = floor (0.5 * NS);
a = 100;
BM = sqrt(2 * a) * exp(1/2);

if strcmp(exc, "bowed")
    FbInit = 1;
else
    FbInit = 0;
end
  
tol = 1e-15;

%% Mass Variables
f1 = 500;    % fundamental frequency [Hz] (< 1 / (k * pi) (< 14,037 Hz))
w1 = 2 * pi * f1;   % angular frequency (< 2 / k (< 88,200 rad/s))
M = 0.001;

u2 = -0.99;
u2Prev = -1;
u2Next = 0;

%% Collision Variables
cL = floor (NS * bridgeLoc) + 1; % bridge location
alpha = 1.3;
K = 5e15;

%% Excitation
amp = -1;

if strcmp(exc, "cos")
    width = 30;
    loc = 1/2;
    startIdx = floor(floor(loc * NS) - width / 2) + 1;
    endIdx = floor(floor(loc * NS) + width / 2) + 1;
    u1Next(startIdx : endIdx) = u1(startIdx : endIdx) + amp * (1 - cos(2 * pi * [0:width]' / width)) / 2;
end
u1 = u1Next;
u1Prev = u1Next;

%% Initialise
psiPrev = 0;

%% Initialise Energy Vectors
kinEnergy1 = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
energy1 = zeros(lengthSound, 1);

energy2 = zeros(lengthSound, 1);

colEnergy1 = zeros(lengthSound, 1);
colEnergy2 = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

vec = 3:NS-2;
eVec = 2:NS-1; % energy vector
Ibr = zeros(NS,1);
Ibr(cL) = 1;
Jbr = Ibr / h;

Ibow = zeros(NS,1);
Ibow(bP) = 1;
Jbow = Ibow / h;

outputPos = floor(NS / 6);
prog = 0;

rampVal = 10000;
eta = u2 - u1(cL);
etaPrev = u2Prev - u1Prev(cL);
etaNext = eta;
for n = 1:lengthSound 
    eta = u2 - u1(cL);
    etaPrev = u2Prev - u1Prev(cL);
        
    %% Update FDSs without connection-force and collision terms 
    strVec = 3:NS-2;
    
    %% NR collision
    b = -2 * eta + etaPrev + k^2 * w1^2 * u2 + T * k^2 / (rho * A * h^2) * (u1(cL+1) - 2 * u1(cL) + u1(cL-1));
    eps = 1;
    ii = 0;

    while eps > tol && ii < 100
        f = subplus(etaNext)^(alpha+1) - subplus(etaPrev)^(alpha+1);
        df = (alpha+1) * subplus(etaNext)^alpha * 0.5 * (1+sign(etaNext));
        gg = etaNext - etaPrev;
        dgg = 1;
        
        coeff = (k^2 / (rho * A * h) + k^2 / M) * K / (alpha+1);
        G = etaNext + coeff * (subplus(etaNext)^(alpha+1)-subplus(etaPrev)^(alpha+1)) / (etaNext - etaPrev) + b;
        Gdiff = 1 + coeff * (df * gg - f * dgg) / gg^2;
        etaNextNew = etaNext - G / Gdiff;
        
        eps = abs(etaNextNew - etaNext);
        etaNext = etaNextNew;
        ii = ii + 1;
    end

    collEffect = K / (alpha + 1) * (subplus(etaNext)^(alpha+1) - subplus(etaPrev)^(alpha+1)) / (etaNext - etaPrev);
    u1Next(strVec) = 2 * u1(strVec) - u1Prev(strVec) + T * k^2 / (rho * A * h^2) * (u1(strVec+1) - 2 * u1(strVec) + u1(strVec-1)) + Jbr(strVec) * k^2 / (rho * A) * collEffect;
    u2Next = 2 * u2 - u2Prev - k^2 * w1^2 * u2 - k^2 / M * collEffect;
    disp(ii + " " + (etaNext - (u2Next - u1Next(cL))))

    %% Calculate energy of the system
    % Energy String
    kinEnergy1(n) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy1(n) = T / 2 * 1/h * sum((u1(3:NS) - u1(2:NS-1)) .* (u1Prev(3:NS) - u1Prev(2:NS-1)));
    energy1(n) = kinEnergy1(n) + potEnergy1(n);
    
    % Energy Mass
    energy2(n) = M / 2 * (1/k * (u2 - u2Prev))^2 + M * w1^2 / 2 * u2 * u2Prev;

    % Collision energy
    colEnergy1(n) = K / (alpha + 1) * 0.5 * (subplus(eta)^(alpha+1) + subplus(etaPrev)^(alpha+1));

%     colEnergy1(n) = K / (alpha + 1) * (subplus(eta)^(alpha+1));
    
    % Total Energy
    totEnergy(n) = energy1(n) + energy2(n) + colEnergy1(n);
        
    %% Update states    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;

    out(n) = u1Next(outputPos);
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
%             legend(["String", "Bridge", "Body"])

            %...plate
            subplot(4,1,2)
            imagesc(reshape(u3Next, [Ny,Nx])')
            xticks([])
            yticks([])
            title("Body (plate)")
            set(gca,'Fontsize', 16)
            subplot(4,1,3);
            
            range = fs / 441;
            if (n > range)
                plot(gSave(n-range:n));
            end
%             ylim([0,50000]);
            title("$g_1^n$", 'interpreter', 'latex')
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
            subplot(3,1,1)
            cla
            hold on;
            %...string
            plot(u1Next, 'Linewidth', 2);    
            %...mass
            scatter(cL, u2Next, 400, '.');    
            % Extra functions
%                 ylim([-1e-4, 1e-4]); % Set y-limit to the amplitude of the raised cosine
            grid on; 
            set(gca, 'Linewidth', 2, 'Fontsize', 16)
            title("State of the system")
%             legend(["String", "Bridge", "Body"])
%             subplot(4,1,2)
%             plot(energy1(1:n));
%             title("String energy",'interpreter', 'latex')
% 
%             subplot(4,1,3)
%             plot(energy2(1:n));
%             title("Mass energy",'interpreter', 'latex')
%             subplot(4,1,4)
%             plot(colEnergy1(1:n));
%             title("Collision energy",'interpreter', 'latex')
            subplot(3,1,2)
            hold off
            plot(energy1(1:n) + energy2(1:n))
            hold on
%             plot(energy2(1:n))
            plot(colEnergy1(1:n))
% plot(totEnergy(1:n))
            subplot(3,1,3)
            plot(totEnergy(1:n) / totEnergy(1) - 1)

%             plot(totEnergy(1:n));
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
plot(g2Save)