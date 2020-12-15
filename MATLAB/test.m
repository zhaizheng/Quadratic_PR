%%
syms x y;
g = [y^2-2*x; -2*y^3+2*x*y];
Tg = [y^2-2*x -2*y^3+2*x*y];
H = [-2 2*y;2*y 2*x-6*y^2];
ans = Tg*g*H*g - Tg*H*g*g;
value = subs(ans,[x,y],[0,1])
for i = 1:2
	t = factor(ans(i))
end
%%

