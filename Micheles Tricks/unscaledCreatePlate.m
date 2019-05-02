function [B, C, N, Nx, Ny, h, kappa, D, Dkappa, DD] = unscaledCreatePlate (Lx, Ly, rho, E, H, s0, s1, k)
    
    Dkappa = E * H^3 / (12 * (1-0.3^2));
    kappa = sqrt(Dkappa / (rho * H));
    LWRatio = Lx/Ly;
    h = 2*sqrt(k*(s1^2+sqrt(kappa^2+s1^2)));

    Nx = floor(Lx/h);
    Ny = floor(Ly/h);
    
    h = max(Lx/Nx, Ly/Ny);
    N = (Nx-1)*(Ny-1);
    
    % generate scheme matrices
    Dxx = sparse(toeplitz([-2;1;zeros(Nx-3,1)]));
    Dyy = sparse(toeplitz([-2;1;zeros(Ny-3,1)]));
    D = kron(speye(Nx-1), Dyy)+kron(Dxx, speye(Ny-1)); 
    DD = D*D/h^4; 
    B = sparse((2 * rho * H / k^2 * speye(N)- Dkappa * DD) / (rho * H / k^2 + s0/k));
    C = -((1-s0*k)/(1+s0*k))*speye(N); % <- not sure if this is still right 
end