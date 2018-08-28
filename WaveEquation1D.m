clear all;
close all;
clc;

fs = 44100;      % Sampling rate
k = 1 / fs;     % Time step
L = 0.4;       % String length
rho = 7850;     % Density of steel [kg/m^3]
c = 200;        % Wave speed
h = c * k;      % Grid spacing

lambdaSq = (c * k / h)^2; % Courant number squared

P = 1/4;              % Plucking position
N = floor(L / h) + 1; % Number of grid-points

%% Raised cosine input
cosWidth = round(N / 10);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);

%% Initialise state vectors
u = zeros(N, 1);
u(floor(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev = u;
uNext = zeros(N, 1);

%% Boundary Conditions
bound = "Dir";
 
%% Extra Settings
lengthSound = fs / 10;
matrix = true;

%% Matrix Representation
if matrix
    if strcmp(bound, "Dir")
        N = N - 2;
        Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
        I = sparse(1:N, 1:N, ones(1, N), N, N);
        N = N + 2;
        range = 2 : N - 1;  
    else
        if strcmp(bound, "Neu")
            Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
                sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
            I = sparse(1:N, 1:N, ones(1, N), N, N);
            range = 1 : N; 
        end
    end

    B = 2 * I + lambdaSq * Dxx;
    out = zeros(lengthSound, 1);
    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
    plotString = false;

    for n = 1 : lengthSound
        uNext(range) = B * u(range) - uPrev(range);
        kinEnergy(n) = 1 / 2 * sum (h * ((1 / k * (u(range) - uPrev(range))).^2));
        potEnergy(n) = c^2 / 2 * sum (1 / h * ...
            (u([1 range] + 1) - u([1 range])) .* (uPrev([1 range] + 1) - uPrev([1 range])));
        out(n) = uNext(PIdx);

        if plotString
            plot(uNext);
            ylim([-1 1]);
            drawnow;
        end

        uPrev = u;
        u = uNext;
    end
    totEnergy = kinEnergy + potEnergy;
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    plot(totEnergy);


else
    %% For-loop representation
    potEnergy = zeros(lengthSound, 1);
    kinEnergy = zeros(lengthSound, 1);

    for n = 1 : lengthSound
        for l = 1 : length(u) - 1
            if l ~= 1
                uNext(l) = 2 * (1 - lambdaSq) * u(l) + lambdaSq * (u(l - 1) + u(l + 1)) - uPrev(l);
            end
                kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
                potEnergy(n) = potEnergy(n) + c^2 / 2 * 1 / h * ...
                    (u(l + 1) - u(l)) * (uPrev(l + 1) - uPrev(l));
        end
        uPrev = u;
        u = uNext;

        out(n) = uNext(PIdx);
    end
    figure(2);
    totEnergy = kinEnergy + potEnergy;
    kinEnergy = (kinEnergy-totEnergy(1))/totEnergy(1);
    potEnergy = (potEnergy-totEnergy(1))/totEnergy(1);
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    % plot(kinEnergy); hold on;
    % plot(potEnergy);
    plot(totEnergy);
end