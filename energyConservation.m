% Draw the paraboloid for different values of a
energy = 200;
for a = -1:0.05:1 
    func = x.^2 + y.^2 + 2 * a * y * x;
    clf;
    mesh(func); hold on;
%     mesh (sqrt(energy/(1-a^2)) * (func < sqrt(energy/(1-a^2))), 'FaceColor','red');

    view(30,60);
    
    drawnow;
end

% Draw slices of the paraboloid from the center outward 
% a = -0.5;
% for ell = 101 :-1: 1 
%     func = x.^2 + y.^2 + 2 * a * y * x;
%     imagesc(round(func) == round(func(ell,ell)));
%     drawnow;
% end