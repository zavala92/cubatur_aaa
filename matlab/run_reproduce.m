function run_reproduce()
%RUN_REPRODUCE Run all numerical experiments and refresh paper figures.
%
%   From the matlab/ directory, run
%       run_reproduce
%
%   This calls run_freeze, which executes every experiment used in the
%   paper, then copies the regenerated figure PDFs into ../paper/figures.
here = fileparts(mfilename('fullpath'));
old = pwd;
cleanup = onCleanup(@() cd(old));
cd(here);

run_freeze();

figs = {'cub2_ellipse.pdf', 'cub3_domains.pdf', 'cub3b_qd.pdf', ...
        'cub4_potential.pdf', 'cub6_perturbed.pdf'};
figdir = fullfile(here, '..', 'paper', 'figures');
if exist(figdir, 'dir') ~= 7
    mkdir(figdir);
end
for k = 1:numel(figs)
    src = fullfile(here, figs{k});
    dst = fullfile(figdir, figs{k});
    if exist(src, 'file') == 2
        copyfile(src, dst);
    else
        warning('run_reproduce:missingFigure', ...
                'Expected figure was not generated: %s', figs{k});
    end
end
fprintf('Reproduction run complete. Figures synced to %s\n', figdir);
end
