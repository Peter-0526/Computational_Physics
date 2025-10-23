tic
N1 = 100000000;
n_total = 3000; % 集团粒子总数
ra = Random_generator_16807(N1);     % 随机数,用于模拟随机游走
rb = Random_generator_16807(200000);     % 随机数,用于放置种子
i_ra = 1;
i_rb = 1;

L = 300;
L0 = 2*L +1; % 网格边长
A = zeros(L0);  % 1:有种子
count = 1;      % 计数器
X_seed = zeros(1000,1); % 模拟种子的运动
Y_seed = zeros(1000,1);

A(L,L) = 1; % 在(L,L)处放置起始点

while count <= n_total
    [Ix,Iy] = find(A);  % 集团坐标值
    d_max = max(sqrt((Ix-L).^2 + (Iy-L).^2));
    r_max = max(d_max+5,20);   % 集团到原点最大距离

    % 种子初始坐标值
    X_seed(1) = round(floor(r_max)*cos(2*pi*rb(i_rb))) + L;
    Y_seed(1) = round(floor(r_max)*sin(2*pi*rb(i_rb))) + L;
    i_rb = i_rb + 1;
    
    for i = 2:N1+1
        if ra(i_ra) < 0.25
            X_seed(i) = X_seed(i-1) + 1;
            Y_seed(i) = Y_seed(i-1);
        elseif  ra(i_ra) < 0.5
            X_seed(i) = X_seed(i-1);
            Y_seed(i) = Y_seed(i-1) + 1;
        elseif  ra(i_ra) < 0.75
            X_seed(i) = X_seed(i-1) - 1;
            Y_seed(i) = Y_seed(i-1);
        else
            X_seed(i) = X_seed(i-1);
            Y_seed(i) = Y_seed(i-1) - 1;
        end
        if (X_seed(i)-L)^2 + (Y_seed(i)-L)^2 > 9*r_max^2      % 走出初始圆半径$\sqrt(3)$倍,舍弃该粒子
            break
        end
        if (X_seed(i) <= 1 || X_seed(i) >= L0 || Y_seed(i) <= 1 || Y_seed(i) >= L0)
            break
        end
        if A(X_seed(i)-1,Y_seed(i)) == 1 || A(X_seed(i)+1,Y_seed(i)) == 1 ...
                || A(X_seed(i),Y_seed(i)-1) == 1 || A(X_seed(i),Y_seed(i)+1) == 1 ...
                || A(X_seed(i)-1,Y_seed(i)-1) == 1 || A(X_seed(i)+1,Y_seed(i)-1) == 1 ...
                || A(X_seed(i)-1,Y_seed(i)+1) == 1 || A(X_seed(i)+1,Y_seed(i)+1) == 1
            A(X_seed(i),Y_seed(i)) = 1;     % 周围有集团粒子,更新A,count
            count = count +1;
            break
        end
        i_ra = i_ra + 1;
    end
end
toc

figure;
h = scatter(Ix,Iy,8,'filled');
axis([0,L0,0,L0])
grid on
grid minor
axis equal
title(['2维DLA模型（粒子总数N = ',num2str(n_total),'）'])