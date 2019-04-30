clear all;
close all;

fs = 44100;
L = 1;
c = 196*2;
k = 1 / fs;     % Time step
h = c * k;      % Grid spacing
BP = 1/pi;      % Bowing Position
lambdaSq = (c * k / h)^2; % Courant number squared
N = floor(L / h) + 1; % Number of grid-points
N = N - 2;
%% Raised cosine input
cosWidth = round(N / 10);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);

%% Initialise state vectors
u = zeros(N, 1);
uPrev = u;
uNext = zeros(N, 1);

u2 = zeros(N, 1);
uPrev2 = u2;
uNext2 = zeros(N, 1);

%% Extra Settings
lengthSound = fs;
drawString = true;

%% Matrix Representation
Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
    sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
    sparse(1:N-1, 2:N, ones(1, N-1), N, N));
I = sparse(1:N, 1:N, ones(1, N), N, N);

uPrev = u;

B = 2 * I + lambdaSq * Dxx;
kappa = 0;
kOh = (kappa/h^2)^2;
gOh = (c/h)^2;
s0 = 0;
s1 = 0;
phi = (2 * s1)/(k * h^2);
bB = (sparse(3:N, 1:N-2, kOh * ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (6*kOh + 2 * gOh - 2/k^2 + 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, kOh * ones(1, N-2), N, N));
bC = sparse(2:N, 1:N-1, phi * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (2/k^2 - 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, phi * ones(1, N-1), N, N);

VbInit = 0.2;
Fb = 50;
a = 100;

BM = sqrt(2 * a) * exp(1/2);
bp = floor(BP * N);
q2 = -VbInit;
% qPrev = -Vb;
qPrev = 0;
qPrev2 = 0;
I = zeros(N,1);
I(bp) = 1;
rampVal = 10000;
out = zeros(lengthSound, 1);
J = 1/h * I;

qPrev = 0;
tol = 1e-7;

for n = 1 : lengthSound
    if n < rampVal
        Vb = VbInit * n / rampVal;
    else
        Vb = VbInit;
    end
    
    % Newton-Raphson
    b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
    eps = 1;
    i = 0;
    while eps>tol
        q=qPrev-(1/h * Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0*qPrev+b)/...
         (1/h * Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0);
        eps = abs(q-qPrev);
        qPrev = q;
        i = i + 1;
        if i > 10000
            disp('Nope');
        end
    end
    qSave(n) = q;
    
%     q = 2 / k * (u(bp) - uPrev(bp)) - qPrev - 2 * Vb;
    %% simple backwards
%     q = (u2(bp) - uPrev2(bp)) / k - Vb;
    
    excitation = k^2 * J * Fb * sqrt(2 * a) * q * exp(-a*q^2 + 1/2);
    uNext = B * u - uPrev - excitation;

    
    q2 = 2 / k * (u2(bp) - uPrev2(bp)) - qPrev2 - 2 * Vb;

    Bm = Fb * (sqrt(2 * a) * exp(-a * q2^2 + 1/2)) / h;
    uNext2 = (B * u2 - uPrev2 + I * Bm * k .* uPrev2 + I * k^2 * Bm * Vb)  ./ (I * k * Bm + 1);

    if mod(n,100) == 0 && drawString == false %&& n > 10000
        subplot(2,1,1);
        plot(u);
        subplot(2,1,2);
        plot(u2);
        drawnow;
    end
    
    qPrev = q;
    uPrev = u;
    u = uNext;
    
    qPrev2 = q2;
    uPrev2 = u2;
    u2 = uNext2;
    out(n) = u(floor(2*N / 3));
    out2(n) = u2(floor(2*N / 3));
end

% totEnergy = kinEnergy + potEnergy;
% totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
% maxTotEnergy = max(totEnergy) - min(totEnergy);
% % plot (totEnergy);
plot(out); 
hold on;
plot(out2)
legend(["Output NR","Output Tricks"])