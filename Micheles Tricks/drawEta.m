figure('Renderer', 'painters', 'Position', [100 100 700 400])
plot(eta,(eta+abs(eta) )/ 2, 'k', 'Linewidth', 2)
grid on
set(gca, 'Fontsize', 16)
ylabel("$[\eta]_+$", 'interpreter', 'latex')
xlabel("$\eta$", 'interpreter', 'latex')
ylim([-0.1,1])