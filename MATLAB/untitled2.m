clear all
format rat
syms x y t a b c d x1
%%%  notation
a = (9 - 42*t + 36*t^2)/(-3*(2*t - 1)*(18*t^2 - 21*t + 4));
b = (6 - 24*t + 18*t^2)/(-3*(2*t - 1)*(18*t^2 - 21*t + 4));
c = 1 /(-3*(2*t - 1)*(18*t^2 - 21*t + 4));
d = (2 - 3*t)/(-3*(2*t - 1)*(18*t^2 - 21*t + 4));

%%%%%%relation between  eigenvalues of m-cell and m+1 cell
s = 36*t^3 - 48*t^2 + 15*t;

a1=(9 - 42*s + 36*s^2)/(-3*(2*s - 1)*(18*s^2 - 21*s + 4));
b1=(6 - 24*s + 18*s^2)/(-3*(2*s - 1)*(18*s^2 - 21*s + 4));
y=0;
%collect(3/4*a + c*(2*x + 2*y*s + (2 - 3*s)*(x + y*s - 3/4 *a1) + 3/4*b1) ,t);

%simple((2 - 3*t)*(3/4 - 3/4 *a) + 3/4*b);
  
collect(-1/4*a + c*(2*x + 2*y*s + (2 - 3*s)*(x + y*s + 1/4 *a1) - 1/4*b1),t);

simple((2 - 3*t)*(-1/4 + 1/4 *a) -1/4*b)
