function cub4_potential()
%CUB4_POTENTIAL Near-singular VOLUME potential over a starfish body:
%       Phi(z0) = int int_Omega f(z) log|z - z0| dA,
%   the 2D Newtonian potential of density f. The dbar-antiderivative
%       W(z, zbar) = (zbar - conj(z0)) (log|z - z0| - 1/2)
%   is single-valued for z0 inside or outside Omega, so one code path
%   evaluates the potential everywhere. Interior targets: nodes cluster
%   AT the target. Exterior targets: nodes cluster at the target's
%   reflection inside the domain (the image point), found automatically.
cfg = cub_config();
t = 2*pi*(0:cfg.M-1)'/cfg.M;
FS = cfg.FS;

rf  = @(s) 1 + 0.3*cos(5*s);
rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s);
zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
zs = zf(t);
f1 = @(z) exp(z);

targets = struct('name', {}, 'z0', {});
tdir = 0.55;                                  % direction for near targets
zb = zf(tdir);                                % boundary point
nrm = 1i*zfp(tdir)/abs(zfp(tdir));            % outward-ish normal
if real(conj(zb)*nrm) < 0, nrm = -nrm; end
targets(1) = struct('name', 'interior, central',      'z0', 0.20 + 0.10i);
targets(2) = struct('name', 'interior, near boundary','z0', zb - 0.05*nrm);
targets(3) = struct('name', 'exterior, delta = 0.05', 'z0', zb + 0.05*nrm);
targets(4) = struct('name', 'exterior, delta = 0.01', 'z0', zb + 0.01*nrm);

results = cell(numel(targets), 1);
for j = 1:numel(targets)
    z0 = targets(j).z0;
    dist = min(abs(zs - z0));
    Wf = @(s) (conj(zf(s)) - conj(z0)).*(log(abs(zf(s) - z0)) - 0.5);
    om = Wf(t);
    Iref = cub_ref(f1, zf, zfp, Wf);
    an = []; ae = []; ee = []; keep = [];
    for deg = 6:6:126
        [zk, wk, epsf] = cub_rule(zs, om, deg, cfg.pdeg);
        if isempty(zk), continue; end
        an(end+1) = numel(zk); %#ok<*AGROW>
        ae(end+1) = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        ee(end+1) = epsf;
        if isempty(keep) && ae(end) < 1e-10, keep = zk; end
    end
    [be, bi] = min(ae);
    fprintf('%-24s dist = %.3f: best rel err %.2e at n = %d (eps %.1e)\n', ...
            targets(j).name, dist, be, an(bi), ee(bi));
    if isempty(keep), keep = zk; end
    results{j} = struct('an', an, 'ae', ae, 'keep', keep, 'z0', z0);
end

figure('Color','w','Position',[100 100 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(real(zs([1:end 1])), imag(zs([1:end 1])), '-', 'Color', cfg.col.curve, ...
     'LineWidth', 1.1); hold on
r2 = results{2}; r3 = results{3};
plot(real(r2.keep), imag(r2.keep), '.', 'Color', cfg.col.node, 'MarkerSize', 11)
plot(real(r2.z0), imag(r2.z0), 'p', 'Color', cfg.col.node, 'MarkerSize', 11, ...
     'MarkerFaceColor', cfg.col.node)
plot(real(r3.keep), imag(r3.keep), '.', 'Color', cfg.col.gauss, 'MarkerSize', 11)
plot(real(r3.z0), imag(r3.z0), 'p', 'Color', cfg.col.gauss, 'MarkerSize', 11, ...
     'MarkerFaceColor', cfg.col.gauss)
axis equal; box on
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
mk = {'o', 's', '^', 'v'};
cl = {cfg.col.node, [0.5 0.2 0.6], cfg.col.gauss, cfg.col.alt};
for j = 1:4
    semilogy(results{j}.an, results{j}.ae, mk{j}, 'Color', cl{j}, ...
             'MarkerSize', 5, 'LineWidth', 0.9); hold on
end
ylim([1e-16 10]); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('number of nodes $n$', 'FontSize', FS)
ylabel('relative error', 'FontSize', FS)
legend({targets.name}, 'FontSize', FS-3, 'Box', 'off', 'Location', 'southwest')
exportgraphics(gcf, fullfile(cfg.here, 'cub4_potential.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end
