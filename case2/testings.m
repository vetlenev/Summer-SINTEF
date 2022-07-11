mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui
A = linspace(20, 50, 200);
A = reshape(A, 20, 10);
%save('summer_sintef/case2/data/test3.mat', 'A');

s = 5;
t = 'Testing';
t1 = 'New %d';

r = sprintf(strcat(t, 'New %d'), s) 