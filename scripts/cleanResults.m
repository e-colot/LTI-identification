function cleanResults(folder, fs, doDelete)
% cleanResults Remove .mat files that do NOT start with 'out[fs]kACQ' or with 'out[fs]kAWGEXC'
% with fs in kHz
%
% Usage:
%   cleanResults()                % dry-run in current folder (lists files to remove)
%   cleanResults(folder)          % dry-run in folder
%   cleanResults(folder,true)     % actually delete the files
%
% Example:
%   cleanResults('C:\data', 5, true)

if nargin < 1 || isempty(folder)
    folder = pwd;
end
if nargin < 2 || isempty(fs)
    fs = 0;
end
if nargin < 3 || isempty(doDelete)
    doDelete = false; % default: dry run
end

if ~isfolder(folder)
    error('Folder not found: %s', folder);
end

files = dir(fullfile(folder, '*.mat'));
files = [files; dir(fullfile(folder, '*.txt'))]; % to remove 'Runner_setup.txt' 
removed = 0;
for k = 1:numel(files)
    name = files(k).name;
    % keep only names that start with 'out[fs]kACQ' or with 'out[fs]kAWGEXC'
    if isempty(regexp(name, ['^out' num2str(fs) 'kACQ.*\.mat$'], 'once'))
        fullpath = fullfile(folder, name);
        fprintf('%s %s\n', ternary(doDelete, 'Deleting:', 'Would delete:'), fullpath);
        if doDelete
            try
                delete(fullpath);
                removed = removed + 1;
            catch ME
                warning('Failed to delete %s: %s', fullpath, ME.message);
            end
        end
    end
end
if doDelete
    fprintf('Deleted %d file(s).\n', removed);
else
    fprintf('Dry run complete. Pass doDelete=true to actually remove files.\n');
end
end

function out = ternary(cond, a, b)
% small helper to avoid inline if-else for printing
if cond
    out = a;
else
    out = b;
end
end