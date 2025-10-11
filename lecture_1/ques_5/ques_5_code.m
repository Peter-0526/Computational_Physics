tic
a1 = 16807;m = 2147483647;q = 127773;r = 2836;
N = 20000000;
I = zeros(1,N);
cl = fix(clock());
I(1) = cl(1) + 70*(cl(2) + 12 * (cl(3) +...
    31 * (cl(4) + 23 * (cl(5) + 59 * cl(6)))));

for i = 2 : N+2
    I(i) = a1 * mod(I(i-1), q) - r * floor(I(i-1)/q);
    if I(i) < 0
        I(i) = I(i) + m;
    end
end
I = I./m;
u = 2.*I(1:N/2)-1;
v = 2.*I(N/2+1:N)-1;

%判断(u,v)点是否满足u^2+v^2 <= 1
r_2 = u.^2 + v.^2;
find_r = find(r_2 <= 1);

X = 2*u(find_r).*sqrt(1-r_2(find_r));
Y = 2*v(find_r).*sqrt(1-r_2(find_r));
Z = 1 - 2.*r_2(find_r);
N_plot = 10000;
scatter3(X(1:N_plot),Y(1:N_plot),Z(1:N_plot),2,'filled','b')
axis equal;
xlabel('X')
ylabel('Y')
zlabel('Z')

figure
h1 = 0.02;
xedg = -1:h1:1;
yedg = -1:h1:1;
h = histogram2(X,Y,xedg,yedg,'Normalization','pdf','FaceColor','flat');
count1 = h.Values;
xlabel('X')
ylabel('Y')
zlabel('p(x,y)')

% 与理论值的对比
figure
[x_t,y_t] = meshgrid(-1+h1/2:h1:1-h1/2,-1+h1/2:h1:1-h1/2);
z_t = (1-x_t.^2-y_t.^2 > 0).*1./(2*pi.*sqrt(1-x_t.^2-y_t.^2));
z_err = (z_t>0).*abs(count1-z_t)./(z_t+eps);
surf(x_t,y_t,z_err,'FaceColor','flat')
ave_err = sum(sum(z_err))/((2/h1)^2.*pi/4);
xlabel('X')
ylabel('Y')
zlabel('error of p(x,y)')
disp(['平均误差值为: ', num2str(ave_err)])

toc