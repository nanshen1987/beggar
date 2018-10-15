% Step 1: arrange points in a matrix format
points = [1 1 1 ; 2 2 2;3 3 3;4 4 4;5 5 5];

% Step 2: find the mean of the points
avg = mean(points, 1);

% Step 3: subtract the mean from all points
subtracted = bsxfun(@minus, points, avg);

% Step 4 : perform SVD
[~, ~, V] = svd(subtracted);

% Step 5: find the direction vector 
%        (which is the right singular vector corresponding to the largest singular value)
direction = V(:, 1);

% Line is 'avg' and 'direction'
p0 = avg;
d = direction;

% Parametric equation: P = p0 + t*d
norm(cross(points(4,:)-p0,d))/norm(d)