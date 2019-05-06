clear all
close all
ratio = zeros(10);
corr = zeros(10);
s0 = 0.001;
s0Save = zeros(10,1);
K1Save = zeros(10,1);
for i = 1:5
    K1 = 100;
    for j = 1:5
        ratio(i, j) = stringSpringMassPlateFunc(s0, K1, 4000);
        if i == 1
            K1Save(j) = K1;
        end
        K1 = K1 * 10;
    end
    s0Save(i) = s0;
    s0 = s0 * 10;
    disp(i);
end
mesh(log(ratio))