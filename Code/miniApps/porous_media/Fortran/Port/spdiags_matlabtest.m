B = zeros(6, 3);
idiags = zeros(1, 3);

% instantiate B
B(1,1) = 1;
B(2,1) = 2;
B(3,1) = 3;
B(4,1) = 4;
B(5,1) = 5;
B(6,1) = 6;

B(1,2) = 7;
B(2,2) = 8;
B(3,2) = 9;
B(4,2) = 10;
B(5,2) = 11;
B(6,2) = 12;

B(1,3) = 13;
B(2,3) = 14;
B(3,3) = 15;
B(4,3) = 16;
B(5,3) = 17;
B(6,3) = 18;

% instantiate idiags
idiags(1) = -2;
idiags(2) = 0;
idiags(3) = 2;

disp('Matrix A');
spdiags(B, idiags, 6, 6)
disp('Matrix B');
spdiags(B, idiags, 6, 5)
disp('Matrix C');
spdiags(B, idiags, 5, 6)
disp('Matrix D');
spdiags(B, idiags, 5, 5)
disp('Matrix E');
spdiags(B, idiags, 5, 4)
disp('Matrix F');
spdiags(B, idiags, 4, 5)
disp('Matrix G');
spdiags(B, idiags, 2, 2)
disp('Matrix H');
spdiags(B, idiags, 1, 1)
