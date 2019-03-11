function z = label_mat2vec(Z)

if size(Z,2) > 1
    [~,z] = max(Z,[],2);
else
    z = Z;
end