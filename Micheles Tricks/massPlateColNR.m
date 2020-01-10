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

%% Mass Variables
f1 = 100;    % fundamental frequency [Hz] (< 1 / (k * pi) (< 14,037 Hz))
omega0 = 2 * pi * f1;   % angular frequency (< 2 / k (< 88,200 rad/s))
M = 0.001;

u2 = offset;
u2Prev = u2 * 1.1;
% u2Prev = -1e-10;
u2Next = 0;

w2 = offset;
w2Prev = w2 * 1.1;
% w2Prev = -1e-10;
w2Next = 0;

%% Plate Variables
Lx = 1.5;
Ly = 0.4;
rhoP = 7850;
EP = 2e11;
H = 0.001;

s0P = 0;
s1P = 0;

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
alpha = 1.3;
K = 5 * 10^4;

%% Excitation

%% Initialise
psiPrevU = 0;

%% Initialise Energy Vectors U
kinEnergy2u = zeros(lengthSound, 1);
potEnergy2u = zeros(lengthSound, 1);
energy2u = zeros(lengthSound, 1);

kinEnergy3u = zeros(lengthSound, 1);
potEnergy3u = zeros(lengthSound, 1);
energy3u = zeros(lengthSound, 1);

colEnergyu = zeros(lengthSound, 1);
totEnergyu = zeros(lengthSound, 1);

%% Initialise Energy Vectors W
kinEnergy2w = zeros(lengthSound, 1);
potEnergy2w = zeros(lengthSound, 1);
energy2w = zeros(lengthSound, 1);

kinEnergy3w = zeros(lengthSound, 1);
potEnergy3w = zeros(lengthSound, 1);
energy3w = zeros(lengthSound, 1);

colEnergyw = zeros(lengthSound, 1);
totEnergyw = zeros(lengthSound, 1);

JbrP = zeros(NP, 1);
JbrP(brP) = 1 / hP^2;
IbrP = zeros(NP, 1);
IbrP(brP) = 1;

prog = 0;

DeltaDelta = Denergy * Denergy;

etawNext = w3(brP) - w2;

colTol = 1e-15;
etaSolut = 0;
for n = 1:lengthSound 
            
    etau = u3(brP) - u2;
    %% Explicitly calculate g
    g = 0;
    if alpha == 1
        if etau > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * subplus(etau)^((alpha - 1)/2);
    end
    gSave(n) = g;
    
    %% Update FDSs without connection-force and collision terms     
    u2Next = (2 * u2 - u2Prev - k^2 * omega0^2 * (u2 - offset));
    u3Next = BP * u3 + CP * u3Prev;
    
    s0P = s0P / (rhoP * H);
    s1P = s1P / (rhoP * H);
    
    %% Matrix solution
    
    Adiv = 1 + g^2*k^2 / (4*M) + g^2*k^2 / (4*hP^2*rhoP*H);
    v = (u3Next(brP) - u3Prev(brP) - u2Next + u2Prev) / (2*k) - k*psiPrevU*g / (2*M) - psiPrevU*g*k / (2*rhoP*H*hP^2);

    solut = v / Adiv;
    
    psiU = g*k*solut + psiPrevU;
    
    u2Next = u2Next + k^2 / M * 0.5 * (psiU + psiPrevU) * g;
    u3Next = u3Next - JbrP * k^2 / (rhoP * H) * 0.5 * (psiU + psiPrevU) * g;

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
    etaw = w3(brP) - w2;
    etawPrev = w3Prev(brP) - w2Prev;
    
    bb = -2 * w3(brP) + w3Prev(brP) + 2 * w2 - w2Prev + D * k^2 / (rhoP * H * hP^4) * DeltaDelta(brP, :) * w3 - k^2 * omega0^2 * (w2 - offset);
    
    eps = 1;
    ii = 0;

    while eps > colTol && ii < 100
        
        f = subplus(etawNext)^(alpha+1) - subplus(etawPrev)^(alpha+1);
        df = (alpha+1) * subplus(etawNext)^alpha * 0.5 * (1 + sign(etawNext));
        gg = etawNext - etawPrev;
        dgg = 1;
        
        coeff = (k^2 / (rhoP * H * hP^2) + k^2 / M) * K / (alpha+1);
        
        G = etawNext + coeff * f / gg + bb;
        dG = 1 + coeff * (df * gg - f * dgg) / gg^2;
        
        etawNextNew = etawNext - G / dG;
        
        eps = abs(etawNextNew - etawNext);
        etawNext = etawNextNew;
        
        ii = ii + 1;
    end
    disp(ii)
    
    force(n) = K / (alpha + 1) * (subplus(etawNext)^(alpha+1) - subplus(etawPrev)^(alpha+1)) / (etawNext - etawPrev);
    
    w2Next = (2 * w2 - w2Prev - k^2 * omega0^2 * (w2 - offset)) + k^2 / M * force(n);
    w3Next = BP * w3 + CP * w3Prev - JbrP * k^2 / (rhoP * H) * force(n);

    %% Calculate energy of the system
    
    % Energy Mass
    kinEnergy2u(n) = M / 2 * (1/k * (u2 - u2Prev)).^2;
    potEnergy2u(n) = M / 2 * omega0^2 * ((u2 - offset) * (u2Prev - offset));
    energy2u(n) = kinEnergy2u(n) + potEnergy2u(n);
    
    % Energy Plate
    kinEnergy3u(n) = ((rhoP * H) / 2) * hP^2 * sum(sum(1/k^2 * (u3 - u3Prev).^2));
    potEnergy3u(n) = D / (2 * hP^2) * sum((Denergy * u3) .* (Denergy * u3Prev));
    energy3u(n) = kinEnergy3u(n) + potEnergy3u(n);
    
    % Connection and Collision energies
    colEnergyu(n) = psiPrevU^2 / 2;
    
    % Total Energy
    totEnergyu(n) = energy2u(n) + energy3u(n) + colEnergyu(n);
    
    %% Calculate energy of the system
    % Energy Mass
    kinEnergy2w(n) = M / 2 * (1/k * (w2 - w2Prev)).^2;
    potEnergy2w(n) = M / 2 * omega0^2 * ((w2 - offset) * (w2Prev - offset));
    energy2w(n) = kinEnergy2w(n) + potEnergy2w(n);
    
    % Energy Plate
    kinEnergy3w(n) = ((rhoP * H) / 2) * hP^2 * sum(sum(1/k^2 * (w3 - w3Prev).^2));
    potEnergy3w(n) = D / (2 * hP^2) * sum((Denergy * w3) .* (Denergy * w3Prev));
    energy3w(n) = kinEnergy3w(n) + potEnergy3w(n);
    
    % Connection and Collision energies
    colEnergyw(n) = K / (alpha + 1) * 0.5 * (subplus(etaw)^(alpha+1) + subplus(etawPrev)^(alpha+1));
    
    % Total Energy
    totEnergyw(n) = energy2w(n) + energy3w(n) + colEnergyw(n);
    
    %% Update states
    psiPrevU = psiU; 

    u2Prev = u2;
    u2 = u2Next;
    
    u3Prev = u3;
    u3 = u3Next;
    
    w2Prev = w2;
    w2 = w2Next;
    
    w3Prev = w3;
    w3 = w3Next;
    
    out2(n) = u3Next(outputPosPlate);
    out3(n) = u2Next;
    diffSave(n) = u2Next - w2Next;
    %% Draw functions
    if mod(n,drawSpeed) == 0 && n >= drawStart && drawThings == true
        subplot(4,2,1)
        hold off;
        scatter(0, u2Next);
        xlim([-1 1])
        ylim([-1e-6 1e-4])
        hold on;
        scatter(0, u3Next(brP))
        
        subplot(4,2,2)
        hold off;
        scatter(0, w2Next);
        xlim([-1 1])
        ylim([-1e-6 1e-4])
        hold on;
        scatter(0, w3Next(brP))
        
        subplot(4,1,2)
        plot(diffSave(1:n))
        
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
        plot(totEnergyw(1:n) / totEnergyw(1) - 1);
%         plot(energy3w(10:n));
%         legend("1", "2", "3", "col1", "col2")
        title("Normalised Energy")

        drawnow;
    end
end
posOut2 = out2 - min(out2);
totOut2 = (posOut2/max(abs(posOut2)) - 0.5) * 2;
posOut3 = out3 - min(out3);
totOut3 = (posOut3/max(abs(posOut3)) - 0.5) * 2;
totOut4 = (g1Save/max(abs(g1Save)) - 0.5) * 2;
totOut = totOut1/5 + totOut2 + totOut3;
% soundsc(totOut1 /5 + totOut2 + totOut3,fs)

plot(gSave)