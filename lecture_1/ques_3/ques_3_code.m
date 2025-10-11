tic

% 产生随机数
a = 16807;m = 2147483647;q = 127773;r = 2836;
N =2000;
I = zeros(1,N);
cl = fix(clock());
I(1) = cl(1) + 70*(cl(2) + 12 * (cl(3) +...
31 * (cl(4) + 23 * (cl(5) + 59 * cl(6)))));

for i = 2 : N+2
    I(i) = a * mod(I(i-1), q) - r * floor(I(i-1)/q);
    if I(i) < 0
       I(i) = I(i) + m;
    end
end
I = I./m;

% 前1000个为第一个坐标值\theta,后1000为\varphi
r = 1;
theta = acos(I(1:N/2));
varphi = 2*pi*I(N/2+1:N);

% 绘制图形
X = r.*sin(theta).*cos(varphi);
Y = r.*sin(theta).*sin(varphi);
Z = r.*cos(theta);

scatter3(X,Y,Z,'filled','o','b');
xlabel('X','FontSize',20);
ylabel('Y','FontSize',20);
zlabel('Z','FontSize',20);
axis equal;
toc