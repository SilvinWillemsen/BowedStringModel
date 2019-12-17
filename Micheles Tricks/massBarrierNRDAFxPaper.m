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

u1Next = zeros(lengthSound, 1);
u1 = -0.999;
u1Prev = -1;

u2Next = zeros(lengthSound, 1);
u2 = -1;
u2Prev = -1;

eta1Prev = u1 - b; 
eta2Prev = u2 - b; 

psiPrev = 0;
bSave = zeros(lengthSound, 1);
bSave = bSave + b;
tol = 1e-10;

eta1NextPrev = u1 - b;
eta1Prev = u1Prev - b;
eta1Next = u1 - b;
u1NextNRPrev = u1 - u1Prev;
for n = 1:lengthSound
    eta1 = u1 - b;
    eta2 = u2 - b;
    
    %% NR
    eps = 1;
    ii = 0;
    
    B = M / k^2 * (-2 * u1 + u1Prev) + M * w0^2 * u1;
    while eps > tol && ii < 100
        f = (subplus(u1NextNRPrev))^(alpha+1) - subplus(u1Prev)^(alpha+1);
        df = (alpha+1) * subplus(u1NextNRPrev)^alpha * 0.5 * (1+sign(u1NextNRPrev));
        gg = u1NextNRPrev - u1Prev;
        dgg = 1;
        G = M/k^2 * u1NextNRPrev + K / (alpha + 1) * (subplus(u1NextNRPrev)^(alpha+1) - subplus(u1Prev)^(alpha+1))/(u1NextNRPrev - u1Prev) + B;
        Gdiff = M/k^2 + K / (alpha+1) * (df * gg - f * dgg) / gg^2;

        u1NextNR = u1NextNRPrev - G / Gdiff;
        eps = abs(u1NextNR - u1NextNRPrev);
        u1NextNRPrev = u1NextNR;
%         end
        ii = ii + 1;
    end
%     ii
%     eta1Next = u1NextNR - b;
    
%     if eta1Next - eta1Prev == 0
%         collEffect = 0;
%     else
%         collEffect = K / (alpha+1) * (mySubplus(eta1Next)^(alpha+1) - mySubplus(eta1Prev)^(alpha+1)) / (eta1Next - eta1Prev);
%     end
%     if collEffect ~= 0 
%         disp("wait");
%     end
    u1Next(n) = u1NextNR;
%     u1Next(n) = (M / k^2 * (2 * u1 - u1Prev) - M * w0^2 * u1 - collEffect) / (M / k^2);

    H1mass(n) = M / 2 * (1/k * (u1 - u1Prev))^2 + M * w0^2 / 2 * u1 * u1Prev;
    H1coll(n) = K / (alpha + 1) * 0.5 * (subplus(u1)^(alpha+1) + subplus(u1Prev)^(alpha+1));
    H1tot(n) = H1mass(n) + H1coll(n);
%     u1Next(n) = (M / k^2 * (2 * u1 - u1Prev) - M * w0^2 * u1 - collEffect) / (M / k^2);
%     eta1Next = u1Next(n) - b;
%     rPrev = u1Next(n) - u1Prev;
    
    %% calculate g
    g = 0;
    if alpha == 1
        if eta2 > 0
            g = sqrt(K * (alpha+1) / 2);
        end
    else
        g = sqrt(K * (alpha+1) / 2) * mySubplus(eta2)^((alpha - 1)/2);
    end
    gSave(n) = g;
    
    u2Next(n) = (M / k^2 * (2 * u2 - u2Prev) - M * w0^2 * u2 + g^2 / 4 * u2Prev - psiPrev * g) / (M / k^2 + g^2/4);
    
    eta2Next = u2Next(n) - b;
    psi = g/2 * (eta2Next - eta2Prev) + psiPrev;
    
    H2mass(n) = M / 2 * (1/k * (u2 - u2Prev))^2 + M * w0^2 / 2 * u2 * u2Prev;
    H2coll(n) = 0.5 * psiPrev^2;
    H2tot(n) = H2mass(n) + H2coll(n);
    
    psiPrev = psi;
    
    u1Prev = u1;
    u1 = u1Next(n);
    
    u2Prev = u2;
    u2 = u2Next(n);
    
    eta1Prev = eta1;
    eta1 = eta1Next;
    eta2Prev = eta2;
    eta2 = eta2Next;

    if drawThings && n > 10
        subplot(2,1,1)
        hold off;
        plot(u1Next(1:n));
        hold on;
        plot(u2Next(1:n));
        plot(bSave(1:n))
        subplot(2,1,2)
%         hold off;
%         plot(H1mass(1:n));
%         hold on;
%         plot(H1coll(1:n));
        plot(H1tot(1:n) / H1tot(1) - 1);
        drawnow;
    end
end