function cub7_details()
%CUB7_DETAILS Algorithmic-detail measurements for the paper.
%   (a) Fit residual vs validation residual: the residual of the
%       least-squares fit is evaluated both on the M fitting samples and
%       on an offset validation grid of 4M points; the ratio measures how
%       much the discrete fit residual underestimates the continuous one.
%   (b) Sensitivity to the polynomial degree L in the fit basis.
%   Domain: the starfish of the paper; area integral (w = 1), f = e^z.
cfg = cub_config();
rf  = @(s) 1 + 0.3*cos(5*s);
rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s);
zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
f1 = @(z) exp(z);
Iref = cub_ref(f1, zf, zfp, @(s) conj(zf(s)));

M = cfg.M;
t  = 2*pi*(0:M-1)'/M;
Mv = 4*M;
tv = 2*pi*((0:Mv-1)' + 0.5)/Mv;               % offset validation grid
zs = zf(t);  om  = conj(zs);
zsv = zf(tv); omv = conj(zsv);

fprintf('(a) fit vs validation residual, L = %d\n', cfg.pdeg);
fprintf('  deg    n    err        eps_fit    eps_val    ratio\n');
worst = 0;
for deg = 30:30:150
    [zk, wk, ef] = cub_rule(zs, om, deg, cfg.pdeg);
    [~,  ~,  ev] = cub_rule(zs, om, deg, cfg.pdeg, zsv, omv);
    e = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
    fprintf('%5d  %3d   %.2e   %.2e   %.2e   %.2f\n', ...
            deg, numel(zk), e, ef, ev, ev/ef);
    worst = max(worst, ev/ef);
end
fprintf('  worst validation/fit ratio: %.2f\n', worst);

fprintf('(b) sensitivity to polynomial degree L (deg = 150)\n');
for L = [4 10 16]
    best = inf;
    for deg = 30:30:150
        [zk, wk, ~] = cub_rule(zs, om, deg, L);
        e = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        best = min(best, e);
    end
    fprintf('  L = %2d:  best rel err %.2e\n', L, best);
end
end
