function cub3b_qd()
%CUB3B_QD Figure for the quadrature-domain experiment: the domain
%   phi(D), phi(zeta) = zeta + 0.3 zeta^2, its exact Aharonov-Shapiro
%   identity pi[(1+2a^2) f(0) + a f'(0)], and the rule discovered by
%   rational approximation, whose confluent node pair realizes the
%   derivative node. Right panel: error against degree; the identity is
%   found at degree 2 and nothing changes afterwards.
cfg = cub_config();
a = 0.3;
phi  = @(zeta) zeta + a*zeta.^2;
phip = @(zeta) 1 + 2*a*zeta;
t = 2*pi*(0:cfg.M-1)'/cfg.M;
zq = phi(exp(1i*t)); omq = conj(zq);
zfun  = @(s) phi(exp(1i*s));
zpfun = @(s) 1i*exp(1i*s).*phip(exp(1i*s));
omfun = @(s) conj(phi(exp(1i*s)));

f1 = @(z) exp(z);
Iref = cub_ref(f1, zfun, zpfun, omfun);
Iexact = pi*((1 + 2*a^2)*f1(0) + a*1);        % f'(0) = 1 for e^z
fprintf('QD: boundary ref vs exact identity: %.2e\n', abs(Iref - Iexact));

dn = []; de = [];
for deg = 2:2:20
    [zk, wk, ~] = cub_rule(zq, omq, deg, cfg.pdeg);
    if isempty(zk), continue; end
    dn(end+1) = deg; %#ok<*AGROW>
    de(end+1) = abs(sum(wk.*f1(zk)) - Iexact)/abs(Iexact);
end
[zk, wk] = cub_rule(zq, omq, 10, cfg.pdeg);
[~, ord] = sort(abs(wk), 'descend');
pair = zk(ord(1:2));
fprintf('confluent pair at %.2e%+.2ei and %.2e%+.2ei (separation %.1e)\n', ...
        real(pair(1)), imag(pair(1)), real(pair(2)), imag(pair(2)), ...
        abs(pair(1) - pair(2)));

FS = cfg.FS;
figure('Color','w','Position',[100 100 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(real(zq([1:end 1])), imag(zq([1:end 1])), '-', 'Color', cfg.col.curve, ...
     'LineWidth', 1.1); hold on
plot(real(pair), imag(pair), '.', 'Color', cfg.col.node, 'MarkerSize', 14)
text(0.06, 0.09, 'confluent pair $\approx$ node of order 2', ...
     'FontSize', FS-2, 'Color', cfg.col.node)
axis equal; box on
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
semilogy(dn, de, 'o-', 'Color', cfg.col.node, 'MarkerSize', 6, ...
         'LineWidth', 1, 'MarkerFaceColor', 'w')
ylim([1e-16 10]); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('AAA degree', 'FontSize', FS)
ylabel('relative error', 'FontSize', FS)
exportgraphics(gcf, fullfile(cfg.here, 'cub3b_qd.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end
