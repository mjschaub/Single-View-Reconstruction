function H = calc_homography(points_one, points_two)

points_one_t = transpose(points_one);
points_two_t = transpose(points_two);

S1 = std(points_one_t);
S2 = std(points_two_t);

T1 = zeros(3, 3);
T1(1, 1) = 1/S1(1);
T1(2, 2) = 1/S1(2);
T1(3, 3) = 1;
points_one_normalized = points_one_t*T1;
points_one_normalized = transpose(points_one_normalized);

T2 = zeros(3, 3);
T2(1, 1) = 1/S2(1);
T2(2, 2) = 1/S2(2);
T2(3, 3) = 1;
points_two_normalized = points_two_t*T2;
points_two_normalized = transpose(points_two_normalized);

A = zeros(4, 9);
for i = 1:4
    A(i*2 - 1, 1) = -points_one_normalized(1, i);
    A(i*2 - 1, 2) = -points_one_normalized(2, i);
    A(i*2 - 1, 3) = -1;
    
    A(i*2 - 1, 7) = points_one_normalized(1, i)*points_two_normalized(1, i);
    A(i*2 - 1, 8) = points_one_normalized(2, i)*points_two_normalized(1, i);
    A(i*2 - 1, 9) = points_two_normalized(1, i);
    
    A(i*2, 4) = -points_one_normalized(1, i);
    A(i*2, 5) = -points_one_normalized(2, i);
    A(i*2, 6) = -1;
    
    A(i*2, 7) = points_one_normalized(1, i)*points_two_normalized(2, i);
    A(i*2, 8) = points_one_normalized(2, i)*points_two_normalized(2, i);
    A(i*2, 9) = points_two_normalized(2, i);    
end

[U, S, V] = svd(A);
h = V(:, end);

H = zeros(3, 3);
for i = 1:3
    H(i, 1) = h(i*3 - 2);
    H(i, 2) = h(i*3 - 1);
    H(i, 3) = h(i*3);
end

H = inv(T2) * H * T1;



