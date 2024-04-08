function [frobeniusSquares, minFroIndices, minFroSquares, lengthsTable] = FrobeniusTable(gb,mis,rotMats,axes,angles)
    % 初始化Frobenius范数的平方数组
    frobeniusSquares = zeros(length(mis), numel(rotMats));
    
    % 初始化用于存储最小Frobenius范数的平方的索引数组
    minFroIndices = zeros(length(mis), 1);
    minFroSquares = zeros(length(mis), 1);
    
    % 创建一个新的数组来保存对应晶界的最小Frobenius范数的平方的轴角对
    alignedAxisAnglePairs = cell(length(mis), 2);
    
    % 初始化一个8x12的表格，用于存储不同轴角对的晶界长度总和
    lengthsTable = zeros(length(axes), length(angles));
    
    % 计算每个晶界与参考旋转矩阵的Frobenius范数的平方，并累加长度
    for i = 1:length(mis)
        actualRotMat = matrix(mis(i));
        frobeniusSquare = zeros(1, numel(rotMats));
        
        for j = 1:numel(rotMats)
            % 计算Frobenius范数的平方
            frobeniusSquare(j) = norm(actualRotMat - rotMats{j}, 'fro')^2;
        end
        
        [minFrobenius, minIndex] = min(frobeniusSquare);
        frobeniusSquares(i,:) = frobeniusSquare;
        minFroIndices(i) = minIndex;
        minFroSquares(i) = minFrobenius;
        
        [axisIndex, angleIndex] = ind2sub([length(axes), length(angles)], minIndex);% 线性索引回下标
        alignedAxisAnglePairs{i, 1} = axes(axisIndex).hkl; % 轴
        alignedAxisAnglePairs{i, 2} = angles(angleIndex); % 角度
    
        % 如果Frobenius范数的平方小于0.01，累加对应轴角对的晶界长度
        if minFrobenius < 0.01
            lengthsTable(axisIndex, angleIndex) = lengthsTable(axisIndex, angleIndex) + gb(i).segLength;
        end
    end
end