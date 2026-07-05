function [zk, wk, epsf, zout] = cub_rule(zs, om, deg, pdeg, zsv, omv)
%CUB_RULE Cubature rule from rational approximation of boundary data.
%
%   The area integral of an analytic density f against a weight w over a
%   domain Omega equals (1/2i) times the contour integral of f * W over
%   gamma = boundary(Omega), where dW/dzbar = w (Cauchy-Pompeiu). Given
%   samples zs of gamma and om = W restricted to gamma, this routine
%   approximates om on gamma by
%       rho(z) = sum_k c_k/(z - z_k)  +  polynomial,
%   using AAA for the pole locations and a least-squares fit for the
%   coefficients. Poles inside Omega become cubature nodes with weights
%   pi*c_k; poles outside Omega and the polynomial integrate to zero
%   against analytic f and serve only to improve the fit. The sup norm
%   of the fit residual on gamma certifies the cubature error:
%       |I - I_n| <= (|gamma|/2) * epsf * max_gamma |f|.
%
%   [zk, wk, epsf, zout] = cub_rule(zs, om, deg, pdeg, zsv, omv)
%     zs   boundary samples (complex column, ordered along gamma)
%     om   values of W on the samples
%     deg  AAA degree
%     pdeg polynomial degree in the least-squares basis
%     zsv, omv  optional VALIDATION grid (finer, offset from zs); if
%          given, epsf is the residual on the validation grid rather
%          than on the fitting samples
%     zk   nodes (poles inside Omega),  wk  weights (= pi * c_k)
%     epsf sup-norm residual (validation grid if provided, else fit grid)
%     zout poles outside Omega (basis helpers, no weight)
if nargin < 5, zsv = []; omv = []; end
ws = warning('off', 'all');
try
    [~, pol] = aaa(om, zs, 'degree', deg, 'sign', 1, 'lawson', 0, ...
                   'tol', 0, 'cleanup', 0);
catch
    [~, pol] = aaa(om, zs, 'degree', deg, 'lawson', 0, ...
                   'tol', 0, 'cleanup', 0);
end
warning(ws);
pol = pol(isfinite(pol));
dmin = min(abs(zs - pol.'), [], 1);
pol = pol(dmin(:) > 1e-8);

in = inpolygon(real(pol), imag(pol), real(zs), imag(zs));
zk = pol(in);
zout = pol(~in);

sc = max(abs(zs));
A = [1./(zs - zk.'), 1./(zs - zout.'), (zs/sc).^(0:pdeg)];
nrm = sqrt(sum(abs(A).^2, 1));
csc = (A./nrm) \ om;
c = csc ./ nrm.';
if isempty(zsv)
    epsf = max(abs(A*c - om));
else
    Av = [1./(zsv - zk.'), 1./(zsv - zout.'), (zsv/sc).^(0:pdeg)];
    epsf = max(abs(Av*c - omv));
end
wk = pi * c(1:numel(zk));
end
