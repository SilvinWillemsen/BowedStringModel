clear all;
close all;

fs = 44100;
k = 1/fs;
b = 0;
lengthSound = fs;

drawThings = true;

%% Mass variables
M = 0.001;
f0 = 500;
w0 = 2 * pi * f0;

%% Collision variables
K = 5e10;
alpha = 1.3;

etaPrev = u2 - b; 
psiPrev = 0;

%% Initialise state variables
u1Next = zeros(lengthSound, 1);
u1 = -0.999;
u1Prev = -1;

u2Next = zeros(lengthSound, 1);
u2 = -1;
u2Prev = -1;

%% Barrier visualisation
bSave = zeros(lengthSound, 1);
bSave = bSave + b;

%% NR tolerance
tol = 1e-13;

NRtype = "DAFx";
if NRtype == "myNR"
    u1NextNRPrev = u1;
else
    u1NextNR = u1;
    r = u1NextNR - u1Prev;
end

for n = 1:lengthSound
    
    %% My NR
    if NRtype == "myNR"
        eps = 1;
        ii = 0;

        B = M / k^2 * (-2 * u1 + u1Prev) + M * w0^2 * u1;
        while eps > tol && ii < 100
            % variables for quotient rule
            f = (subplus(u1NextNRPrev))^(alpha+1) - subplus(u1Prev)^(alpha+1);
            df = (alpha+1) * subplus(u1NextNRPrev)^alpha * 0.5 * (1+sign(u1NextNRPrev));
            gg = u1NextNRPrev - u1Prev;
            dgg = 1;
            
            G = M/k^2 * u1NextNRPrev + K / (alpha + 1) * (subplus(u1NextNRPrev)^(alpha+1) - subplus(u1Prev)^(alpha+1))/(u1NextNRPrev - u1Prev) + B;
            Gdiff = M/k^2 + K / (alpha+1) * (df * gg - f * dgg) / gg^2;

            % NR for u1Next
            u1NextNR = u1NextNRPrev - G / Gdiff;
            eps = abs(u1NextNR - u1NextNRPrev);
            u1NextNRPrev = u1NextNR;
            ii = ii + 1;
        end
        disp("NR iterations: " + ii)
        u1Next(n) = u1NextNR;
    else
        
    %% DAFx NR
        eps = 1;
        ii = 0;

        a = b - u1Prev;
        B = -2 * u1 + 2 * u1Prev + k^2 * w0^2 * u1;
        while eps > tol && ii < 100
            % variables for quotient rule
            f = (mySubplus(r-a))^(alpha+1) - mySubplus(-a)^(alpha+1);
            df = (alpha+1) * mySubplus(r-a)^alpha * 0.5 * (1+sign(r-a));
            gg = r;
            dgg = 1;
            
            G = r + k^2/M * K/(alpha+1) * (f / gg) + B;
            Gdiff = 1 + k^2 / M * K / (alpha+1) * (df * gg - f * dgg) / gg^2;

            % NR for r (u1Next - u1Prev)
            rNew = r - G / Gdiff;
            u1NextNR = rNew + u1Prev;
            eps = abs(rNew - r);
            r = rNew;
            ii = ii + 1;
        end
        disp("NR iterations: " + ii)
        u1Next(n) = u1NextNR;
    end
    
    % energy
    H1mass(n) = M / 2 * (1/k * (u1 - u1Prev))^2 + M * w0^2 / 2 * u1 * u1Prev;
    H1coll(n) = K / (alpha + 1) * 0.5 * (subplus(u1)^(alpha+1) + subplus(u1Prev)^(alpha+1));
    H1tot(n) = H1mass(n) + H1coll(n);
    
    %% Non-iterative methods
    
    % calculate g
    eta2 = u2 - b;
    g = 0;
    if alpha == 1
        if eta2 > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * mySubplus(eta2)^((alpha - 1)/2);
    end
    gSave(n) = g;
    
    % calculate system
    u2Next(n) = (M / k^2 * (2 * u2 - u2Prev) - M * w0^2 * u2 + g^2 / 4 * u2Prev - psiPrev * g) / (M / k^2 + g^2/4);
    
    eta2Next = u2Next(n) - b;
    psi = g/2 * (eta2Next - etaPrev) + psiPrev;
    
    % energy
    H2mass(n) = M / 2 * (1/k * (u2 - u2Prev))^2 + M * w0^2 / 2 * u2 * u2Prev;
    H2coll(n) = 0.5 * psiPrev^2;
    H2tot(n) = H2mass(n) + H2coll(n);
    
    %% Update states
    psiPrev = psi;
    
    u1Prev = u1;
    u1 = u1Next(n);
    
    u2Prev = u2;
    u2 = u2Next(n);
    
    etaPrev = eta2;
    eta2 = eta2Next;

    %% visualise stuff
    if drawThings && n > 10
        subplot(3,1,1)
        hold off;
        plot(u1Next(1:n));
        hold on;
        plot(u2Next(1:n));
        plot(bSave(1:n));
        legend("NR", "Non-it", "Barrier");
        title("State of the system")
        subplot(3,1,2)
        plot(H1tot(1:n) / H1tot(1) - 1);
        title("Normalised energy NR")
        subplot(3,1,3)
        plot(H2tot(1:n) / H2tot(1) - 1);
        title("Normalised energy Non-iterative method")
        drawnow;
    end
end