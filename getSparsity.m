function S = getSparsity(totalSize)
    e = ones(totalSize,1);
    A = spdiags([e e e],-1:1,totalSize,totalSize);
    S = full(A);
end