clear all;
close all;
clc
maxTotEnergyVect = zeros(1000,1);
BC = ["clamped" "ss" "free"];
fs = zeros(1000,1);
L = zeros(1000,1);
c = zeros(1000,1);

for i = 1:1000
    bound(i) = BC(randi([1 3]));
    fs(i) = round(4000 + rand*40000);
    L(i) = 0.1 + 0.9 * rand;
    c(i) = round(20 + rand*500);
    maxTotEnergyVect(i) = IdealBarEquation (fs(i), L(i), bound(i)); 
end
plot(maxTotEnergyVect);