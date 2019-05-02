function [B, C, N, Nx, Ny, h, kappa, D, Denergy] = newCreatePlate (Lx, Ly, rho, H, E, s0, s1, k)
    D = E * H^3 / (12 * (1 - 0.3^2));
    kappa = sqrt(D / (rho * H));
    LWRatio = Lx/Ly;
    h = 2*sqrt(k*(s1^2+sqrt(kappa^2+s1^2)));

    Nx = floor(Lx/h);
    Ny = floor(Ly/h);
    
    h = max(Lx/Nx, Ly/Ny);
    N = (Nx-1)*(Ny-1);
    
    % generate scheme matrices
    Dxx = sparse(toeplitz([-2;1;zeros(Nx-3,1)]));
    Dyy = sparse(toeplitz([-2;1;zeros(Ny-3,1)]));
    Denergy = kron(speye(Nx-1), Dyy)+kron(Dxx, speye(Ny-1)); 
    DD = Denergy*Denergy/h^4; 
    B = sparse((2*speye(N)-kappa^2*k^2*DD)/(1+s0*k));
    C = -((1-s0*k)/(1+s0*k))*speye(N); 
    
end