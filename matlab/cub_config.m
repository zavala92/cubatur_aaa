function cfg = cub_config()
%CUB_CONFIG Shared configuration for the cubature-from-rational-
%   approximation project (2D generalization; standalone project).
here = fileparts(mfilename('fullpath'));
cfg.here = here;
add_dependency_paths(here);
require_dependency_functions();
cfg.deps.aaa = which('aaa');
cfg.deps.jacpts = which('jacpts');

cfg.M = 800;            % boundary samples fed to AAA
cfg.pdeg = 10;          % polynomial columns in the least-squares fit

cfg.FS = 14;
cfg.col.node  = [0.85 0.10 0.10];
cfg.col.curve = [0.00 0.35 0.72];
cfg.col.gauss = [0.15 0.15 0.15];
cfg.col.alt   = [0.00 0.35 0.72];
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 13);
set(groot, 'defaultAxesLineWidth', 0.6);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesGridAlpha', 0.08);
end

function add_dependency_paths(here)
% Optional local hook. Copy local_paths_template.m to local_paths.m and
% edit it for machine-specific paths if automatic discovery is not enough.
hook = fullfile(here, 'local_paths.m');
if exist(hook, 'file') == 2
    run(hook);
end

if exist('jacpts', 'file') ~= 2 || exist('aaa', 'file') ~= 2
    add_existing_paths(chebfun_candidates(here));
end

if exist('aaa', 'file') ~= 2
    add_existing_paths(aaa_candidates(here));
end
end

function require_dependency_functions()
missing = {};
if exist('aaa', 'file') ~= 2
    missing{end+1} = 'aaa'; %#ok<AGROW>
end
if exist('jacpts', 'file') ~= 2
    missing{end+1} = 'jacpts'; %#ok<AGROW>
end
if ~isempty(missing)
    error(['Missing required MATLAB functions: %s.\n' ...
           'Install Chebfun and add it to the MATLAB path, set ' ...
           'CHEBFUN_PATH/AAA_PATH, or create local_paths.m from ' ...
           'local_paths_template.m.'], strjoin(missing, ', '));
end
end

function paths = chebfun_candidates(here)
home = user_home();
paths = {getenv('CHEBFUN_PATH'), getenv('CHEBFUN_ROOT'), ...
         fullfile(here, 'chebfun'), ...
         fullfile(here, 'extern', 'chebfun'), ...
         fullfile(here, '..', 'chebfun')};
if ~isempty(home)
    paths = [paths, { ...
        fullfile(home, 'Documents', 'MATLAB', 'chebfun'), ...
        fullfile(home, 'Documents', 'MATLAB', 'surfacefun', 'chebfun'), ...
        fullfile(home, 'MATLAB', 'chebfun')}]; %#ok<AGROW>
end
end

function paths = aaa_candidates(here)
home = user_home();
paths = {getenv('AAA_PATH'), getenv('AAA_ROOT'), ...
         fullfile(here, 'aaa'), ...
         fullfile(here, 'extern', 'aaa')};
if ~isempty(home)
    paths = [paths, { ...
        fullfile(home, 'MATLAB', 'aaa'), ...
        fullfile(home, 'MATLAB_local', 'MATLAB', 'aaa')}]; %#ok<AGROW>
end
end

function add_existing_paths(paths)
for k = 1:numel(paths)
    p = paths{k};
    if ~isempty(p) && exist(p, 'dir') == 7
        addpath(p);
    end
end
end

function home = user_home()
home = getenv('HOME');
if isempty(home)
    home = getenv('USERPROFILE');
end
end
