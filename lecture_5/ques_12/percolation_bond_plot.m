init_plot
figure

t = tiledlayout(3,5);   % 绘图布局结构
n = 4;

bonds = zeros(6,1);
for i = 1:n
   bonds(i) = 1; 
end
x1 = [[1,1];[1,2];[1,3];[2,1];[2,2];[2,3]];
z1 = [[2,3];[5,6];[2,5];[1,2];[4,5];[1,4]];
c = 1;

while c <= 720/(factorial(n)*factorial(6-n))
    next_tile(x1)
    for i = 1:6
        if bonds(i) == 1
            plot(x1(z1(i,:),1),x1(z1(i,:),2),'r')
        else
            plot(x1(z1(i,:),1),x1(z1(i,:),2),'k')
        end
    end
    c = c+1;
    bonds = next_bond(bonds,n);
end
title(t,['b = 2,通路n = ',num2str(n),'方阵键逾渗示意图'],'FontSize',18);

function init_plot
    [x,y] = meshgrid(1:20,1:20);
    scatter(x,y,10,'b','filled')
    % 网格线
    hold on
    grid on
    ax = gca;
    xticks(0:20)
    yticks(0:20)
    ax.GridLineWidth = 1;
    ax.GridColor = 'k';
    x1 = 2.5:2:18.5;

    % 分块，重整化群
    for i = 1:length(x1)
        plot([0,20],[x1(i),x1(i)],'--','Color','r')
        plot([x1(i),x1(i)],[0,20],'--','Color','r')
    end
    xlabel('X','FontSize',18)
    ylabel('Y','FontSize',18)
    title('键逾渗重整化群示意图','FontSize',18)
end

function next_tile(x1)
    nexttile
    hold on
    scatter(x1(:,1),x1(:,2),20,'blue','filled')
    axis([0,3,0,4])
    xticks(0:3)
    yticks(0:3)
    grid on
end

function b = next_bond(bonds,n)

    l = length(bonds);
    if n == 0 && n == l     % 全通路/全不通
        b = bonds;
        return
    end    
    c = l;        
    b = bonds;
    if bonds(l) == 0        % bond(end) == 0,将最后一个1后移一位
        while c >= 1 && bonds(c) == 0   
            c = c-1;
        end
        b(c) = 0;
        b(c+1) = 1;
        return
    else                    % bond(end) == 1,去掉所有后缀1,将剩余最后一个1后移一位
        while bonds(c) == 1
            c = c-1;
        end
        n_1 = l - c;
        while c >= 1 && bonds(c) == 0
            c = c-1;
        end
        if c >= 1           % 有n-1个后缀1,第一个1后移一位,其余1接在后面
            if n_1 == n-1
                b = zeros(6,1);
                for i = 1:n
                    b(c+i) = 1;
                end
                return
            else            % 仅改变去掉后缀后最后一个1,后缀的1接上
                b(c) = 0;
                for i = 1:n_1+1
                    b(c+i) = 1;
                end
                for i = n_1+2:l-c
                    b(c+i) = 0;
                end
                return
            end
        else
            b = zeros(6,1);
            for i = 1:n
                b(i) = 1; 
            end
            return
        end
    end
end

