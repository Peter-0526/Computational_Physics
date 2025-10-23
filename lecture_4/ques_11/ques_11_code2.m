N1 = 5000;
X_ra = Random_generator_16807(N1);  % 随机数
i_ra = 1;
N_max = 20000; % 微分方程最多迭代次数
n_total = 100; % 集团粒子总数
L = 200;
L0 = 2*L +3; % 网格边长
A = zeros(L0);
count = 1;      % 计数器
tol = 1e-5;     % Laplace方程tolerance
w = 1.9;

A(L,L) = 1; % 在(L,L)处放置起始点
A0 = A;     % 数值为1:对应点放置种子
eta = 2;

tic
while count < n_total
    for k = 1:N_max
        A1 = A;
        for i= 3:2*L+1
            for j = 3:2*L+1
                if A0(i,j) == 0
                    A(i,j) = A(i,j)*(1-w) + w*(A(i-1,j) + A(i+1,j) + A(i,j-1) + A(i,j+1))/4;
                end
            end
        end
        if abs(max(max(A - A1))) < tol
            break;
        end
    end
    if k >= N_max
        disp('发散')
        break;
    end
    % 集团边缘随机选择一个位置产生新粒子
    [A_addx,A_addy] = find((A0(3:2*L+1,3:2*L+1)-1).*(A0(2:2*L,2:2*L) + A0(3:2*L+1,2:2*L) + A0(4:2*L+2,2:2*L) ...
        + A0(4:2*L+2,3:2*L+1) + A0(4:2*L+2,4:2*L+2) + A0(3:2*L+1,4:2*L+2) ...
        + A0(2:2*L,4:2*L+2) + A0(2:2*L,3:2*L+1)));
    A_addx = A_addx+2;
    A_addy = A_addy+2;
    P = zeros(length(A_addx),1);
    for j = 1:length(A_addx)
        P(j) = (1 - A(A_addx(j),A_addy(j)))^eta;
    end
    X_add = find_index(P,X_ra(count));
    A0(A_addx(X_add),A_addy(X_add)) = 1;
    A(A_addx(X_add),A_addy(X_add)) = 1;
    count = count+1;
end

toc

[Ix,Iy] = find(A0);
figure;
h = scatter(Ix,Iy,8,'filled');
axis([0,L0,0,L0])
grid on
grid minor
axis equal
title(['介电击穿模型（粒子总数N = ',num2str(n_total),',\eta = ',num2str(eta),'）'])