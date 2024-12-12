function import_dir(dir_name, option)

    % Set default value for dir_name if not provided
    if nargin < 1
        dir_name = 'utils';
    end

    % Set default value for option if not provided
    if nargin < 2
        option = 'load';
    end

    if strcmp(option, 'load')
        % Define the folder path relative to the current directory
        folderPath = fullfile(pwd, dir_name);

        % Check if the folder is already on the path
        if ~contains(path, folderPath)
            % Add the folder to the path
            addpath(folderPath);
        end
    elseif strcmp(option, 'unload')
        if isfolder(dir_name)
            % Remove the directory from the path
            rmpath(dir_name);
        else
            disp(['Directory "', dir_name, '" is not in the MATLAB search path.']);
        end
    else
        disp(['Error: "', option, '" is not a valid option.']);
    end

end
