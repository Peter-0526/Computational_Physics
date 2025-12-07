N1 = 100000000;
n_total = 5000; % 集团粒子总数
ra = Random_generator_16807(N1);     % 随机数,用于模拟随机游走
rb = Random_generator_16807(200000);     % 随机数,用于放置种子
i_ra = 1;
i_rb = 1;
L = 300;

[L0, Ix, Iy] = DLA_2D(L, n_total, rb, i_rb, N1, ra, i_ra);

DLA_plot(Ix, Iy, L0, n_total);
figure;

% SandBox法
R = sqrt((Ix-L).^2 + (Iy-L).^2);
R_max = max(R);
k = 20; 
X1 = linspace(log(R_max/k),0.8*log(R_max),k)';
N1 = zeros(k,1);
for i = 1:k
    N1(i) = length(find(R <= exp(X1(i))));
end
LX1 = [ones(k,1),X1];
Ln1 = log(N1);
b = LX1\Ln1;
% 拟合精度
LnCalc1 = LX1*b;
Rsq = 1 - sum((Ln1 - LnCalc1).^2)/sum((Ln1 - mean(Ln1)).^2);
disp(['D = ',num2str(b(2)),' Rsq = ',num2str(Rsq)])
plot(X1,log(N1),'o-')
xlabel('$\log(r)$','Interpreter', 'latex','FontSize',18,'FontName', 'Arial')
ylabel('$\log(N)$','Interpreter', 'latex','FontSize',18,'FontName', 'Arial')
title('Sandbox','Interpreter', 'latex','FontSize',18,'FontName', 'Arial')
figure;


% 密度-密度相关函数法
Delta = min([max(Ix-L),max(L-Ix),max(Iy-L),max(L-Iy)]);
X2 = log(round(exp((linspace(log(Delta/k),0.9*log(Delta),k)'))));
N2 = zeros(k,1);
for m = 1:k
    for i = 1:5000
        flag = zeros(1,4);
        for j = 1:5000
            if sum(flag) == 4
                break;
            end
            if flag(1)==0 && Ix(i)-round(exp(X2(m))) == Ix(j) && Iy(i) == Iy(j)
            N2(m) = N2(m) + 1/4;
            flag(1) = 1;
            end
            if flag(2)==0 && Ix(i)+round(exp(X2(m))) == Ix(j) && Iy(i) == Iy(j)
            N2(m) = N2(m) + 1/4;
            flag(2) = 1;
            end
            if flag(3)==0 && Iy(i)-round(exp(X2(m))) == Iy(j) && Ix(i) == Ix(j)
            N2(m) = N2(m) + 1/4;
            flag(3) = 1;
            end
            if flag(4)==0 && Iy(i)+round(exp(X2(m))) == Iy(j) && Ix(i) == Ix(j)
            N2(m) = N2(m) + 1/4;
            flag(4) = 1;
            end
        end
    end
end
LX2 = [ones(k,1),X2];
Ln2 = log(N2);
b2 = LX2\Ln2;
% 拟合精度
LnCalc2 = LX2*b2;
Rsq2 = 1 - sum((Ln2 - LnCalc2).^2)/sum((Ln2 - mean(Ln2)).^2);
disp(['D = ',num2str(2+b2(2)),' Rsq = ',num2str(Rsq2)])
plot(X2,log(N2),'o-')
xlabel('$\log(r)$','Interpreter', 'latex','FontSize',18,'FontName', 'Arial')
ylabel('$\log(\rho)$','Interpreter', 'latex','FontSize',18,'FontName', 'Arial')
title('密度-密度相关函数法','FontSize',18)


function [L0, Ix, Iy] = DLA_2D(L, n_total, rb, i_rb, N1, ra, i_ra)  % 生成集团粒子坐标，中心点(L,L)
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
        if (X_seed(i)-L)^2 + (Y_seed(i)-L)^2 > 3*r_max^2      % 走出初始圆半径\sqrt(3)倍,舍弃该粒子
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
end

function DLA_plot(Ix, Iy, L0, n_total) % 绘制集团粒子
figure;
h = scatter(Ix,Iy,8,'filled');
axis([0,L0,0,L0])
grid on
grid minor
axis equal
title(['2维DLA模型（粒子总数N = ',num2str(n_total),'）'])
end