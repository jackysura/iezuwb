function norms = norm_vectors(vectors)
len = length(vectors);
norms = zeros(len,1);
for i = 1:len
    norms(i) = norm(vectors(i,:));
end