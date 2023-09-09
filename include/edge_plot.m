function contour_plot = edge_plot(A)
A = round(A);
if size(A,1)>size(A,2)
    transpose = true;
    A = A';
else
    transpose = false;
end

for r = sort(unique(A(1,:)))
    if exist('contour_plot','var')
        contour_plot = [contour_plot, [r,r;min(A(2,A(1,:)==r)),max(A(2,A(1,:)==r))]];
    else
        contour_plot = [r,r;min(A(2,A(1,:)==r)),max(A(2,A(1,:)==r))];
    end
end

for c = sort(unique(A(2,:)))
    contour_plot = [contour_plot, [min(A(1,A(2,:)==c)),max(A(1,A(2,:)==c));c,c]];
end

contour_plot = unique(contour_plot','rows');
if ~transpose
    contour_plot = contour_plot';
end