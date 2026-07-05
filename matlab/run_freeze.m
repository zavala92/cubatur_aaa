function run_freeze()
%RUN_FREEZE Preprint-freeze run: every experiment in the paper, in
%   order, in one fresh session. Console output is the record.
fprintf('=== FREEZE RUN %s ===\n', datestr(now));
fprintf('MATLAB %s\n', version);
fns = {@cub1_disk, @cub2_ellipse, @cub3_domains, @cub3b_qd, ...
       @cub3c_resolution, @cub4_potential, @cub4b_resolution, ...
       @cub5_kernels, @cub6_perturbed, @cub7_details};
for k = 1:numel(fns)
    nm = func2str(fns{k});
    fprintf('\n========== %s ==========\n', nm);
    t0 = tic;
    fns{k}();
    fprintf('---- %s done in %.1f s\n', nm, toc(t0));
    close all
end
fprintf('=== FREEZE COMPLETE ===\n');
end
