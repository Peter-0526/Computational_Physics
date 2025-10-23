function idx = find_index(A, x)
idx = zeros(size(x));
A = A/sum(A);
C = cumsum(A);
for i = 1:length(x)
    l = length(C);
    left = 1;  % MATLAB 数组索引从1开始
    right = l;

    while right > left
        mid = floor((left + right - 1) / 2);  % 使用floor确保整数索引

        if x(i) == C(mid)
            idx(i) = mid + 1;
            return;
        elseif x(i) < C(mid)
            right = mid;
        else
            left = mid + 1;
        end
    end

    % 如果循环结束还没找到，返回最后的mid值
    mid = floor((left + right - 1) / 2);
    idx(i) = mid + 1;
end