function remove_subfolders(parent)
%REMOVE_SUBFOLDERS Remove all subfolders of a given folder from MATLAB path,
% but keep the parent folder itself on the path.
%
% Usage:
%   remove_subfolders('/path/to/replication_package')

    % Get all folders including parent and all subfolders
    allPaths = genpath(parent);
    
    % Split into a cell array
    allPathsCell = strsplit(allPaths, pathsep);
    
    % Remove the first entry (the parent folder itself)
    subFolders = allPathsCell(2:end);

    % Remove only subfolders that are actually on the path
    for i = 1:numel(subFolders)
        if ~isempty(strfind(path, subFolders{i}))
            rmpath(subFolders{i});
        end
    end
end

