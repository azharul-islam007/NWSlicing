function [vec_o] = cell2vec_f(cell_i)
% A utilities function to transform cell matrix to a column vector
vec_len = 0;
nr = numel(cell_i);
for r = 1 : nr
    vec_len = vec_len + numel(cell_i{r});
end

vec_o = zeros(vec_len, 1);
next_count = 1;
for r = 1 : nr
    for c = 1 : numel(cell_i{r})
        vec_o(next_count) = cell_i{r}{c};
        next_count = next_count + 1;
    end
end

end