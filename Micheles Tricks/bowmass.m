% matlab script bowmass.m
% finite difference scheme for a bowed mass-spring system
% soft friction characteristic w/iterative Newton-Raphson method

clear all;
close all;
clc;

SR = 44100; % sample rate (Hz)
TF = 10; % simulation duration (s)
drawThings = true;
drawSpeed = 100;
drawStart = 0;%SR * TF / 2;
viewWindow = drawSpeed * 100;
drawFunc = 0;
%%%%%% begin global parameters

K = 1000;
M = 0.001;
R = 5;
w0 = sqrt(K/M);
FBInit = 5; % bow force
vBInit = 0.1; % bow velocity (m/s)
sig = 100; % friction law free parameter (1/m^2) 
tol = 1e-4; % tolerance for Newton-Raphson method
%%%%%% end global parameters


% derived parameters
NF = floor(TF*SR);
k = 1/SR;
A = exp(1/2)*sqrt(2*sig);

% initialize time series/iterative method

u = zeros(NF,1);
u2 = zeros(NF,1);
u3 = zeros(NF,1);
f = zeros(NF,1);
vr = zeros(NF,1);
qlast = 0;
qlast2 = 0;
% % time step restrictions
% if(k>min(1/(pi*f0),exp(1)/(FB*sqrt(2*sig))))
%     error("Time step too large");
% end
%%%%%% start main loop
ramp = 0;

if ramp == 0
    q2Prev = -vBInit;
    q3Prev = -vBInit;
else
    q2Prev = 0;
    q3Prev = 0;
end

for n=3:NF
    if n < ramp
        vB = n * vBInit / ramp;
        FB = n * FBInit / ramp;
    else
        vB = vBInit;
        FB = FBInit;
    end
    % Newton-Raphson method to determine relative velocity
%     b = -(2 * M/k^2) * (u(n-1)-u(n-2)) + (2 * M / k) * vB + K * u(n-1);
    b = 2 * M / k * (vB - 1/k * (u(n-1) - u(n-2))) + K * u(n-1) + R * vB;
    eps = 1;
    i = 0;
    while eps>tol
        q=qlast-(FB*A*qlast*exp(-sig*qlast^2) + (2*M/k + R) * qlast+b)/...
         (FB*A*(1-2*sig*qlast^2)*exp(-sig*qlast^2)+2*M/k+R);
        eps = abs(q-qlast);
        qlast = q;
        i = i + 1;
        if i > 10000
            disp('whut')
        end
    end
    qSave(n) = q;
    u(n) = (M / k^2 * (2*u(n-1) - u(n-2)) - K * u(n-1) + R / (2*k) * u(n-2) - FB * sqrt(2*sig) * q * exp(-sig*q^2 + 1/2)) / (M / k^2 + R / (2*k));
    
    shouldBe0(n) = q - ((u(n) - u(n-2)) / (2*k) - vB);
    
    rOCkinEnergy(n) = M/(2*k^3) * (u(n)-u(n-2)) * (u(n)-2*u(n-1)+u(n-2));
    rOCpotEnergy(n) = (K * u(n-1)) * (u(n) - u(n-2)) / (2*k);
    rOCdampEnergy(n) = R / (4*k^2) * (u(n) - u(n-2))^2;
    rOCbowEnergy(n) = FB * (u(n) - u(n-2)) / (2*k) * sqrt(2*sig) * q * exp(-sig*q^2 + 1/2);
    rOCtotEnergy(n) = rOCkinEnergy(n) + rOCpotEnergy(n) + rOCbowEnergy(n) + rOCdampEnergy(n);
    
    if drawThings && n > viewWindow && drawFunc == 1 && mod(n,drawSpeed) == 0 && n > drawStart
        subplot(3,1,1)
        plot(u(1:n));
        subplot(3,1,2)
        plot(rOCkinEnergy(n-viewWindow:n));
        hold on;
        plot(rOCpotEnergy(n-viewWindow:n));
        plot(rOCbowEnergy(n-viewWindow:n));
        plot(rOCtotEnergy(n-viewWindow:n));
        hold off;
        subplot(3,1,3)
        plot(rOCtotEnergy(n-viewWindow:n));
        drawnow;
    end
    
    q2 = 2 / k * (u2(n-1) - u2(n-2)) - q2Prev - 2 * vB;

    % update position of mass and relative bow velocity
    
    BM = FB * sqrt(2*sig) * exp(-sig * q2^2+1/2);
    u2(n) = (M / k^2 * (2 * u2(n-1) - u2(n-2)) - K * u2(n-1) + R / (2*k) * u2(n-2) ...
        + BM / (2*k) * u2(n-2) + BM * vB) / (M / k^2 + R / (2*k) + BM / (2 * k));
    shouldBe0_2(n) = q2 - ((u2(n) - u2(n-2)) / (2*k) - vB);
    
    rOCkinEnergy2(n) = M/(2*k^3) * (u2(n)-u2(n-2)) * (u2(n)-2*u2(n-1)+u2(n-2));
    rOCpotEnergy2(n) = -(K * u2(n-1)) * (u2(n) - u2(n-2)) / (2*k);
    rOCdampEnergy2(n) = - R / (4*k^2) * (u2(n) - u2(n-2))^2;
    rOCbowEnergy2(n) = -FB * ((u2(n) - u2(n-2)) / (2*k) - vB) * sqrt(2 * sig) * exp(-sig * q2^2 + 1/2) *  (u2(n) - u2(n-2)) / (2*k);
    rOCtotEnergy2(n) = rOCkinEnergy2(n) - rOCpotEnergy2(n) - rOCbowEnergy2(n) - rOCdampEnergy2(n);
%     subplot(2,1,1)
%     plot(u2(1:n))
%     subplot(2,1,2)
    
    if drawThings && n > viewWindow && drawFunc == 2 && mod(n,drawSpeed) == 0 && n > drawStart
        subplot(3,1,1)
        plot(u2(1:n));
        subplot(3,1,2)
        plot(rOCkinEnergy2(n-viewWindow:n));
        hold on;
        plot(rOCpotEnergy2(n-viewWindow:n));
        plot(rOCbowEnergy2(n-viewWindow:n));
        plot(rOCbowEnergy2(n-viewWindow:n));
        hold off;
        subplot(3,1,3)
        plot(rOCtotEnergy2(n-viewWindow:n));
        drawnow;
    end
%     
    if drawThings && n > viewWindow && drawFunc == 0 && mod(n,drawSpeed) == 0 && n > drawStart
%         subplot(4,2,1)
%         plot(u(1:n));
%         title("State u_1")
%         subplot(4,2,2)
%         plot(u2(1:n));
%         title("State u_2")
%         subplot(4,2,3)
%         plot(rOCkinEnergy(n-viewWindow:n));
%         hold on;
%         plot(rOCpotEnergy(n-viewWindow:n));
%         plot(rOCbowEnergy(n-viewWindow:n));
%         hold off;
%         title("Energy")
%         subplot(4,2,4)
%         plot(rOCkinEnergy2(n-viewWindow:n));
%         hold on;
%         plot(rOCpotEnergy2(n-viewWindow:n));
%         plot(rOCbowEnergy2(n-viewWindow:n));
%         hold off;
%         title("Energy")
%         subplot(4,2,5)
%         plot(rOCtotEnergy(n-viewWindow:n));
%         title("Rate-of-change Energy")
%         subplot(4,2,6)
%         plot(rOCtotEnergy2(n-viewWindow:n));
%         title("Rate-of-change Energy")
%         subplot(4,2,7)
%         plot(shouldBe0)
%         subplot(4,2,8)
%         plot(shouldBe0_2)
%         if n > 10000
%             cla
%             plot(u((n-10000):n))
%             hold on
%             plot(u2(n-10000:n))
%             drawnow;
%         end
    end
%     testFunc2 = k^2 * FB * sqrt(2*sig) * q3 * exp(-sig * q3^2+1/2);
%     u3(n) = 2 * u3(n-1) - u3(n-2) - k^2*w0^2 * u3(n-1) - testFunc2;
%     
    q2Prev = q2;
%     q3Prev = q3;
    q2Save(n) = q2;
%     q3Save(n) = q3;
end
%%%%%% end main loop

% plot mass displacement and relative bow velocity
% tax = [0:NF-1]*k;
% subplot(2,1,1); plot(tax, u, "k"); title("Displacement of Mass");
% xlabel("time (s)"); 
% subplot(2,1,2); 
% plot(tax, vr, "k");
% title("Relative Bow Velocity");
% xlabel("time (s)");
timeVecUnderSamp = (0:length(u) - 1) / SR;
% timeVec = (0:length(u) - 1) / SR;
subplot(3,1,1)
plot(timeVecUnderSamp, u); hold on; plot(timeVecUnderSamp, u2); plot(timeVecUnderSamp, u3)
% load saveTest.mat
plot(timeVec, saveNR)
legend(["NR","TRICKS","Backwards", "OversampledNR"])
subplot(3,1,[2,3])
range = 2600:2750;
plot(timeVecUnderSamp(range), u(range)); 
hold on; 
plot(timeVecUnderSamp(range), u2(range)); 
plot(timeVecUnderSamp(range), u3(range));
plot(timeVec(range*10), saveNR(range*10))
freqAxOS = ((0:length(saveNR) - 1) * SR * 10) / length(saveNR);

freqAx = ((0:length(u) - 1) * SR) / length(u);

legend(["NR","TRICKS","Backwards", "OversampledNR"])