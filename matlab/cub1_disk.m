function cub1_disk()
%CUB1_DISK Sanity demo: on a disk the boundary data zbar = R^2/z + c-terms
%   is rational, so the schema should return a SINGLE node at the center
%   with weight pi*R^2 -- the mean value property, discovered numerically.
cfg = cub_config();
t = 2*pi*(0:cfg.M-1)'/cfg.M;

fprintf('--- unit disk ---\n');
zs = exp(1i*t);
[zk, wk, epsf] = cub_rule(zs, conj(zs), 12, cfg.pdeg);
[~, km] = max(abs(wk));
fprintf('nodes: %d;  dominant node = %.2e%+.2ei;  weight - pi = %.2e;  eps = %.1e\n', ...
        numel(zk), real(zk(km)), imag(zk(km)), abs(wk(km) - pi), epsf);
fs = {@(z) exp(z), @(z) 1./(3 - z), @(z) cos(2*z)};
for k = 1:3
    In = sum(wk .* fs{k}(zk));
    Iref = pi * fs{k}(0);                       % mean value property
    fprintf('  f%d: error = %.2e\n', k, abs(In - Iref));
end

fprintf('--- disk of radius 0.7 centered at 0.5 + 0.2i ---\n');
z0 = 0.5 + 0.2i; R = 0.7;
zs = z0 + R*exp(1i*t);
[zk, wk, epsf] = cub_rule(zs, conj(zs), 12, cfg.pdeg);
[~, kmain] = max(abs(wk));
fprintf('nodes: %d;  main node - center = %.2e;  weight - pi R^2 = %.2e;  eps = %.1e\n', ...
        numel(zk), abs(zk(kmain) - z0), abs(wk(kmain) - pi*R^2), epsf);
In = sum(wk .* exp(zk));
fprintf('  exp: error vs pi R^2 exp(center) = %.2e\n', abs(In - pi*R^2*exp(z0)));
end
