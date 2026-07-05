function I = cub_ref(fh, zfun, zpfun, omfun)
%CUB_REF Reference value of the area integral via the boundary identity,
%   I = (1/2i) * contour integral of f * W over gamma, computed with the
%   periodic trapezoidal rule (spectrally accurate for smooth data) at
%   two resolutions, which must agree to 1e-12.
I1 = trap(2^13); I2 = trap(2^14);
assert(abs(I1 - I2) < 1e-12 * max(1, abs(I2)), ...
       'reference self-check failed: %.2e', abs(I1 - I2));
I = I2;
    function I = trap(M)
        t = 2*pi*(0:M-1)'/M;
        I = (1/(2i)) * (2*pi/M) * sum(fh(zfun(t)) .* omfun(t) .* zpfun(t));
    end
end
