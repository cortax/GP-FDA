function component = resample_components( component, Z, x, y, verbose )
    parfor k = 1:size(Z,2)
        sub_y_idx = find(Z(:,k));
        component{k}.ytr = y(:, sub_y_idx);
        component{k}.N = length(sub_y_idx);
        component{k} = nsgpnuts(component{k}, 0.02, 3, verbose);
    end
end
