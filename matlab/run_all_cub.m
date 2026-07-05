function run_all_cub()
%RUN_ALL_CUB Driver for the cubature-from-rational-approximation demos.
fns = {@cub1_disk, @cub2_ellipse, @cub3_domains, @cub4_potential};
for k = 1:numel(fns)
    nm = func2str(fns{k});
    fprintf('\n========== %s ==========\n', nm);
    t0 = tic;
    try
        fns{k}();
        fprintf('---- %s done in %.1f s\n', nm, toc(t0));
    catch err
        fprintf(2, '---- %s FAILED: %s\n', nm, err.message);
        for f = 1:min(3, numel(err.stack))
            fprintf(2, '     at %s line %d\n', err.stack(f).name, ...
                    err.stack(f).line);
        end
    end
    close all
end
fprintf('CUBALLDONE\n');
end
