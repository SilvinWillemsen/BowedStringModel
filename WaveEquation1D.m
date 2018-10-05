function [maxTotEnergy] = WaveEquation1D (fs, L, c, BC)
% fs = 15600;
% L = 1;
% c = 100;
% BC = "Dir";
k = 1 / fs;     % Time step
h = c * k;      % Grid spacing

lambdaSq = (c * k / h)^2; % Courant number squared

% P = 1/4;              % Plucking position
N = floor(L / h) + 1; % Number of grid-points

%% Raised cosine input
% cosWidth = round(N / 10);
% raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);

%% Initialise state vectors
% u = zeros(N, 1);
% u(floor(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
u = rand(N,1);
% uPrev = rand(N,1);
uPrev = u;
uNext = zeros(N, 1);

%% Extra Settings
lengthSound = fs;
drawString = false;

%% Matrix Representation
if strcmp(BC, "Dir")
    N = N - 2;
    Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
        sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
        sparse(1:N-1, 2:N, ones(1, N-1), N, N));
    I = sparse(1:N, 1:N, ones(1, N), N, N);
    N = N + 2;
    u(1) = 0;
    u(N) = 0;
    uPrev = u;
end
if strcmp(BC, "Neu1")
    Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 1], N, N) + ...
                sparse(1:N, 1:N, [-1 -2 * ones(1, N-2) -1], N, N) + ...
                sparse(1:N-1, 2:N, [1 ones(1, N-2)], N, N));
    I = sparse(1:N, 1:N, ones(1, N), N, N);
end
if strcmp(BC, "Neu2")
    Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
                sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
    I = sparse(1:N, 1:N, ones(1, N), N, N);
    scalefac = ones(N,1);
    scalefac(1) = 0.5;
    scalefac(N) = 0.5;
end

B = 2 * I + lambdaSq * Dxx;

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

if strcmp(BC, "Dir")
    for n = 1 : lengthSound
        uNext(2:N-1) = B * u(2:N-1) - uPrev(2:N-1);
        kinEnergy(n) = 1 / 2 * sum (h * ((1 / k * (u - uPrev)).^2));
        potEnergy(n) = c^2 / 2 * sum (1 / h * ...
            (u(2:N) - u(1:N-1)) .* (uPrev(2:N) - uPrev(1:N-1)));
        
        if mod(n,2) == 0 && drawString
            plot(u);
            drawnow;
        end
        
        uPrev = u;
        u = uNext;
    end
end
if strcmp(BC, "Neu1")
    for n = 1 : lengthSound
        uNext = B * u - uPrev;
        kinEnergy(n) = 1 / 2 * sum (h * ((1 / k * (u - uPrev)).^2));
        potEnergy(n) = c^2 / 2 * sum (1 / h * ...
            (u(2:N) - u(1:N-1)) .* (uPrev(2:N) - uPrev(1:N-1)));
        
        if mod(n,2) == 0 %&& drawString
            plot(u)
            drawnow;
        end
        
        uPrev = u;
        u = uNext;
    end
end

if strcmp(BC, "Neu2")
    for n = 1 : lengthSound
        uNext = B * u - uPrev;
        kinEnergy(n) = 1 / 2 * sum (h * scalefac .* ((1 / k * (u - uPrev)).^2));
        potEnergy(n) = c^2 / 2 * sum (1 / h * ...
            (u(2:N) - u(1:N-1)) .* (uPrev(2:N) - uPrev(1:N-1)));
        
        if mod(n,1) == 0 && drawString
%             totEnergy = kinEnergy + potEnergy;
%             totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
%             plot(totEnergy(1:n))
            plot(u);
            drawnow;
        end
        
        uPrev = u;
        u = uNext;
    end
end

totEnergy = kinEnergy + potEnergy;
totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
maxTotEnergy = max(totEnergy) - min(totEnergy);
% plot (totEnergy);
