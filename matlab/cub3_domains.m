function cub3_domains()
%CUB3_DOMAINS General domains. (A) A starfish domain: the nodes trace the
%   analytic skeleton (mother body) of the domain, and the rule converges
%   geometrically for analytic densities; cross-checked against a 2D
%   adaptive integral in polar coordinates. (B) A quadrature domain
%   phi(D), phi(zeta) = zeta + 0.3 zeta^2: the Schwarz function is
%   meromorphic, an EXACT finite quadrature identity exists
%   (Aharonov--Shapiro), and the schema finds it.
cfg = cub_config();
t = 2*pi*(0:cfg.M-1)'/cfg.M;
FS = cfg.FS;

% ---------- (A) starfish ----------
rf  = @(s) 1 + 0.3*cos(5*s);
rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s);
zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
zs = zf(t); om = conj(zs);

f1 = @(z) exp(z);
f2 = @(z) 1./(1.6 - z);
I1 = cub_ref(f1, zf, zfp, @(s) conj(zf(s)));
I2 = cub_ref(f2, zf, zfp, @(s) conj(zf(s)));
Iarea_re = integral2(@(s,r) real(f1(r.*exp(1i*s))).*r, 0, 2*pi, ...
                     0, @(s) rf(s), 'AbsTol', 1e-11, 'RelTol', 1e-11);
Iarea_im = integral2(@(s,r) imag(f1(r.*exp(1i*s))).*r, 0, 2*pi, ...
                     0, @(s) rf(s), 'AbsTol', 1e-11, 'RelTol', 1e-11);
Iarea = Iarea_re + 1i*Iarea_im;
fprintf(['starfish: boundary ref vs integral2: %.2e ' ...
         '(real %.2e, imag %.2e)\n'], ...
        abs(I1 - Iarea), abs(real(I1) - Iarea_re), ...
        abs(imag(I1) - Iarea_im));

an = []; e1 = []; e2 = []; keep = [];
for deg = 6:6:150
    [zk, wk, ~] = cub_rule(zs, om, deg, cfg.pdeg);
    if isempty(zk), continue; end
    an(end+1) = numel(zk); %#ok<*AGROW>
    e1(end+1) = abs(sum(wk.*f1(zk)) - I1)/abs(I1);
    e2(end+1) = abs(sum(wk.*f2(zk)) - I2)/abs(I2);
    if an(end) >= 60 && isempty(keep), keep = zk; end
end
[be, bi] = min(e1);
fprintf('starfish: best rel err %.2e at n = %d\n', be, an(bi));

% ---------- (B) quadrature domain ----------
phi  = @(zeta) zeta + 0.3*zeta.^2;
phip = @(zeta) 1 + 0.6*zeta;
zq  = phi(exp(1i*t));
omq = conj(zq);
Iq1 = cub_ref(f1, @(s) phi(exp(1i*s)), ...
              @(s) 1i*exp(1i*s).*phip(exp(1i*s)), @(s) conj(phi(exp(1i*s))));
[zkq, wkq, epsq] = cub_rule(zq, omq, 10, cfg.pdeg);
[~, ord] = sort(abs(wkq), 'descend');
zkq = zkq(ord); wkq = wkq(ord);
big = abs(wkq) > 1e-10;
fprintf('quadrature domain: %d nodes carry weight (eps = %.1e)\n', sum(big), epsq);
for k = find(big).'
    fprintf('  node %.15f%+.15fi   weight %.15f%+.2ei\n', ...
            real(zkq(k)), imag(zkq(k)), real(wkq(k)), imag(wkq(k)));
end
fprintf('  sum of weights - area(pi(1+2a^2)) = %.2e\n', ...
        abs(sum(wkq) - pi*(1 + 2*0.3^2)));
maxe = 0;
fs = {f1, f2, @(z) z.^3, @(z) z.^7, @(z) cos(3*z)};
for k = 1:numel(fs)
    Ik = cub_ref(fs{k}, @(s) phi(exp(1i*s)), ...
         @(s) 1i*exp(1i*s).*phip(exp(1i*s)), @(s) conj(phi(exp(1i*s))));
    maxe = max(maxe, abs(sum(wkq.*fs{k}(zkq)) - Ik)/max(1, abs(Ik)));
end
fprintf('  exactness over 5 analytic test densities: max rel err = %.2e\n', maxe);
a = 0.3;
fprintf('  confluent-pair check: sum(w*z) - pi*a = %.2e  (identity pi[(1+2a^2)f(0) + a f''(0)])\n', ...
        abs(sum(wkq.*zkq) - pi*a));

% ---------- figure ----------
figure('Color','w','Position',[100 100 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(real(zs([1:end 1])), imag(zs([1:end 1])), '-', 'Color', cfg.col.curve, ...
     'LineWidth', 1.1); hold on
plot(real(keep), imag(keep), '.', 'Color', cfg.col.node, 'MarkerSize', 11)
axis equal; box on
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
semilogy(an, e1, 'o', 'Color', cfg.col.node, 'MarkerSize', 5.5, 'LineWidth', 1); hold on
semilogy(an, e2, 's', 'Color', cfg.col.alt, 'MarkerSize', 5.5)
ylim([1e-16 10]); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('number of nodes $n$', 'FontSize', FS)
ylabel('relative error', 'FontSize', FS)
text(60, 3e-3, '$f = e^z$', 'Color', cfg.col.node, 'FontSize', FS-1)
text(20, 1e-10, '$f = 1/(1.6-z)$', 'Color', cfg.col.alt, 'FontSize', FS-1)
exportgraphics(gcf, fullfile(cfg.here, 'cub3_domains.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end
