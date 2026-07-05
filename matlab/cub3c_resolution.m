function cub3c_resolution()
%CUB3C_RESOLUTION Starfish area integral: floor vs boundary sampling.
cfg = cub_config();
rf  = @(s) 1 + 0.3*cos(5*s);
rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s);
zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
f1 = @(z) exp(z);
Iref = cub_ref(f1, zf, zfp, @(s) conj(zf(s)));
for M = [800 1600 3200]
    t = 2*pi*(0:M-1)'/M;
    zs = zf(t); om = conj(zs);
    best = inf; bn = 0;
    for deg = 10:10:150
        [zk, wk, ~] = cub_rule(zs, om, deg, cfg.pdeg);
        if isempty(zk), continue; end
        e = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        if e < best, best = e; bn = numel(zk); end
    end
    fprintf('M = %4d:  best rel err %.2e at n = %d\n', M, best, bn);
end
end
