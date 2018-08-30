fs = 44100;
c = 100;
L = 1;
k = 1 / fs;     % Time step
h = c * k;      % Grid spacing

lambdaSq = (c * k / h)^2; % Courant number squared

P = 1/4;              % Plucking position
N = floor(L / h) + 1; % Number of grid-points

%% Raised cosine input
cosWidth = round(N / 10);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);

%% Initialise state vectors
u = zeros(N, 1);
u(floor(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev = u;
uNext = zeros(N, 1);

%% Boundary Conditions
bound = "Neu";
 
%% Extra Settings
lengthSound = fs / 10;

%% Matrix Representation
matrix = true;
if matrix
    if strcmp(bound, "Dir")
        Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
        I = sparse(1:N, 1:N, ones(1, N), N, N); 
    else
        if strcmp(bound, "Neu")
            Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
                sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
            Dxxxx = Dxx * Dxx;
            I = sparse(1:N, 1:N, ones(1, N), N, N);
        end
    end

    B = 2 * I + lambdaSq * Dxx;
    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);

    for n = 1 : lengthSound
        uNext(2:N-1) = B(2:N-1, 2:N-1) * u(2:N-1) - uPrev(2:N-1);
        kinEnergy(n) = 1 / 2 * sum (h * ((1 / k * (u(1:N) - uPrev(1:N))).^2));
        potEnergy(n) = c^2 / 2 * sum (1 / h * ...
            (u(2:N) - u(1:N-1)) .* (uPrev(2:N) - uPrev(1:N-1)));
        uPrev = u;
        u = uNext;
    end
    totEnergy = kinEnergy + potEnergy;
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    maxTotEnergy = max(totEnergy) - min(totEnergy);
    plot (totEnergy);
    
else
    i = 1;
    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
    totEnergy = zeros(lengthSound,1);
    for n = 1 : lengthSound
        for l = 1 : length(u)
            if strcmp(bound, "Dir")
                if l ~= N
                    if l ~= 1
                        uNext(l) = 2 * (1 - lambdaSq) * u(l) + lambdaSq * (u(l - 1) + u(l + 1)) - uPrev(l);
                    end
                    potEnergy(n) = potEnergy(n) + c^2 / 2 * 1 / h * ...
                        (u(l + 1) - u(l)) * (uPrev(l + 1) - uPrev(l));
                end
            else
                if strcmp(bound, "Neu")
                    if l == N
                        uNext(N) = 2 * u(N) - uPrev(N) + lambdaSq * (-u(N)+u(N-1));
%                         potEnergy(n) = potEnergy(n) + c^2 / 2 * 1 / h * ...
%                             (u(N - 1) - u(N - 2)) * (uPrev(N - 1) - uPrev(N - 2));
                    else
                        if l == 1
                            uNext(1) = 2 * u(1) - uPrev(1) + lambdaSq*(u(2)-u(1));
                        else
                            uNext(l) = 2 * u(l) + lambdaSq * (u(l - 1) - 2 * u(l) + u(l + 1)) - uPrev(l);
                        end
%                             + (h / (2 * k)) * (u(1) - uPrev(1)) + (h / (2 * k)) * (u(N) - uPrev(N))
                        potEnergy(n) = potEnergy(n) + c^2 / 2 * 1 / h * ...
                                 (u(l + 1) - u(l)) * (uPrev(l + 1) - uPrev(l));
                    end
                    
                end
            end
            kinEnergy(n) = kinEnergy(n) + 1 / 2 * (h / k^2) * (u(l) - uPrev(l))^2;
%             potEnergy(n) = potEnergy(n) + c^2 / 2 * (h * (1 / 2 * (1 / h * ...
%                     (u(l+1) - u(l) - uPrev(l+1) + uPrev(l))))^2 - ...
%                     k^2 / 4 * (h * (1/k * (1/h * (u(l+1) - u(l) - uPrev(l+1) + uPrev(l))))^2));

        end
        
        subplot(2,1,1);
        plot(uNext);
        totEnergy = kinEnergy + potEnergy;
        subplot(2,1,2);
        plot((totEnergy(1:n)-totEnergy(1))/totEnergy(1));
        drawnow;
        uPrev = u;
        u = uNext;
    end
    % totEnergy = kinEnergy + potEnergy;
    totEnergy = kinEnergy + potEnergy;
    % totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    % plot(totEnergy);
    % hold on;
    plot(totEnergy);
    %     hold on; plot(potEnergy2);
    %     hold on; plot(kinEnergy2);
    % plot(out);
end
