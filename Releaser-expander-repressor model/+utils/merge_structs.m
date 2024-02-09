function merged = merge_structs(s1, s2)
    % merge two structures; overwrite duplicates with s2's value
    merged = s1;
    f = fieldnames(s2);
    for i = 1 : numel(f)
        merged.(f{i}) = s2.(f{i});
    end
end
