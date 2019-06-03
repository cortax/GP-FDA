function x = tocolumn(x)
    if ~isvector(x)
        error('Input must be a vector');
    end
    if ~iscolumn(x)
        x = x';
    end
end

