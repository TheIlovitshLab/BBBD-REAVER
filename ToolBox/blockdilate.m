function dilated = blockdilate(A,se)
% Custom version of imdilate that utilizes the fact that only a small
% fraction of the image A contains non-black pixels.
[row,col] = find(A);
min_row = max([min(row)-size(se.Neighborhood,1),1]);
min_col = max([min(col)-size(se.Neighborhood,2),1]);
max_row = min([max(row)+size(se.Neighborhood,1),size(A,1)]);
max_col = min([max(col)+size(se.Neighborhood,2),size(A,2)]);
dil = imdilate(A(min_row:max_row,min_col:max_col),se);
dilated = false(size(A));
dilated(min_row:max_row,min_col:max_col) = dil;
end