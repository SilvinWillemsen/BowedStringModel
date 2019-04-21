plot(u1, 'Linewidth', 2);
hold on;
scatter(cL, u2, 400, '.');
ylim([-0.15 0.15])
grid on
plot([cL-1, cL+1], [b b], 'Linewidth', 5);
set(gca,'Linewidth', 2, 'Fontsize', 16)
legend(["String", "Bridge", "Body"])