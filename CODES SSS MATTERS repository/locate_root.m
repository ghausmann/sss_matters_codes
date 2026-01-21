% locate_root.m
function root = locate_root(start_dir)
    dir_now = start_dir;
    while true
        if exist(fullfile(dir_now,'root_marker.txt'),'file')  % or any file that marks root
            root = dir_now;
            return;
        end
        parent = fileparts(dir_now);
        if strcmp(parent, dir_now)
            error('Cannot locate replication package root.');
        end
        dir_now = parent;
    end
end
