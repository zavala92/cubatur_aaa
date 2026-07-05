function cub6_perturbed()
%CUB6_PERTURBED The classical rule needs the exact geometry; the AAA rule
%   does not. Domain: an ellipse (a, b) = (1, 0.6) perturbed by
%   p*exp(3it). Gauss--Gegenbauer on the focal segment of the underlying
%   ellipse (the classical rule one would use) floors at O(p); the AAA
%   cubature rule, built from the true boundary data, reaches machine
%   precision regardless of p.
cfg = cub_config();
a = 1; b = 0.6; c = sqrt(a^2 - b^2);
t = 2*pi*(0:cfg.M-1)'/cfg.M;
f1 = @(z) exp(z);

ps = [1e-2, 1e-4];
res = cell(numel(ps), 1);
for j = 1:numel(ps)
    p = ps(j);
    zf  = @(s) a*cos(s) + 1i*b*sin(s) + p*exp(3i*s);
    zfp = @(s) -a*sin(s) + 1i*b*cos(s) + 3i*p*exp(3i*s);
    zs = zf(t); om = conj(zs);
    Iref = cub_ref(f1, zf, zfp, @(s) conj(zf(s)));

    gn = 1:24; ge = zeros(size(gn));
    for n = gn
        [xg, wg] = jacpts(n, 0.5, 0.5);
        ge(n) = abs(2*a*b*(wg*f1(c*xg)) - Iref)/abs(Iref);
    end
    an = []; ae = [];
    for deg = 2:2:44
        [zk, wk, ~] = cub_rule(zs, om, deg, cfg.pdeg);
        if isempty(zk), continue; end
        e = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        if isfinite(e), an(end+1) = numel(zk); ae(end+1) = e; end %#ok<*AGROW>
    end
    fprintf('p = %.0e: Gauss-Gegenbauer floor %.2e; AAA best %.2e at n = %d\n', ...
            p, min(ge), min(ae), an(find(ae == min(ae), 1)));
    res{j} = struct('p', p, 'gn', gn, 'ge', ge, 'an', an, 'ae', ae);
end

FS = cfg.FS;
figure('Color','w','Position',[100 100 1000 380]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:2
    r = res{j};
    nexttile
    semilogy(r.gn, r.ge, 's', 'Color', cfg.col.gauss, 'MarkerSize', 5.5); hold on
    semilogy(r.an, r.ae, 'o', 'Color', cfg.col.node, 'MarkerSize', 5.5, ...
             'LineWidth', 1)
    yline(r.p, ':', 'Color', [.5 .5 .5], 'LineWidth', 1);
    ylim([1e-17 10]); xlim([0 44]); grid on; box on
    set(gca, 'YMinorGrid', 'off', 'FontSize', FS-1)
    xlabel('number of nodes $n$', 'FontSize', FS)
    if j == 1
        ylabel('relative error', 'FontSize', FS)
        text(14, 8e-2, 'Gauss--Gegenbauer', 'Color', cfg.col.gauss, ...
             'FontSize', FS-2)
        text(23, 1e-11, 'AAA cubature', 'Color', cfg.col.node, ...
             'FontSize', FS-2)
    end
end
exportgraphics(gcf, fullfile(cfg.here, 'cub6_perturbed.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end
