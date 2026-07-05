function cub4b_resolution()
%CUB4B_RESOLUTION Confirm that the delta = 0.01 floor in cub4 is set by
%   the boundary sampling density, not by the method: rerun the exterior
%   delta = 0.01 target with M = 3200 boundary samples.
cfg = cub_config();
rf  = @(s) 1 + 0.3*cos(5*s);
rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s);
zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
tdir = 0.55;
zb = zf(tdir);
nrm = 1i*zfp(tdir)/abs(zfp(tdir));
if real(conj(zb)*nrm) < 0, nrm = -nrm; end
z0 = zb + 0.01*nrm;
f1 = @(z) exp(z);
Wf = @(s) (conj(zf(s)) - conj(z0)).*(log(abs(zf(s) - z0)) - 0.5);
Iref = cub_ref(f1, zf, zfp, Wf);
for M = [800 1600 3200]
    t = 2*pi*(0:M-1)'/M;
    zs = zf(t); om = Wf(t);
    best = inf; bn = 0;
    for deg = 20:10:150
        [zk, wk, ~] = cub_rule(zs, om, deg, cfg.pdeg);
        if isempty(zk), continue; end
        e = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        if e < best, best = e; bn = numel(zk); end
    end
    fprintf('M = %4d:  best rel err %.2e at n = %d\n', M, best, bn);
end
end
