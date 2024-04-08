function [NormQuatDiff, minQuatnIndices, minQuatnAngles, lengthsTable] = NormQuaTable(gb, mis, axes, angles, axisAnglePairs)
    % 计算长度
    Length = length(axes) * length(angles);
    numMis = length(mis);

    % 初始化矩阵和数组
    NormQuatDiff = zeros(numMis, Length);
    minQuatnIndices = zeros(numMis, 1);
    minQuatnAngles = zeros(numMis, 1);
    lengthsTable = zeros(length(axes), length(angles));

    % 提取mis中的四元数并转换为矩阵
    quatMis = zeros(numMis, 4);
    for i = 1:numMis
        quatMis(i, :) = [mis(i).a, mis(i).b, mis(i).c, mis(i).d];
    end
    
    % 提取所有轴角对的四元数
    quatAxisAngles = zeros(Length, 4);
    for j = 1:Length
        aA = axisAnglePairs{j, 3}; % 假设axisAnglePairs的每个元素都是一个结构，其中包含四元数
        quatAxisAngles(j, :) = [aA.a, aA.b, aA.c, aA.d];
    end

    % 计算所有晶界与每个轴角对之间的距离
    for i = 1:numMis
        actualQuat = repmat(quatMis(i, :), Length, 1); % 重复当前四元数以匹配quatAxisAngles的尺寸
        diffs = actualQuat - quatAxisAngles; % 计算差值
        normQuatDiff = sqrt(sum(diffs.^2, 2))'; % 计算欧几里得距离（二范数）
        
        [minQuatNorm, minQuatnIndice] = min(normQuatDiff);
        NormQuatDiff(i, :) = normQuatDiff;
        minQuatnIndices(i) = minQuatnIndice;
        minQuatnAngles(i) = minQuatNorm;

        [axisIndex, angleIndex] = ind2sub([length(axes), length(angles)], minQuatnIndice);
        % 假设axes是一个结构体数组，angles是一个向量
        % 这里不再保存轴角对，因为它在简化示例中未被使用

        % 更新lengthsTable
        if minQuatNorm < 0.05 % 你可能需要定义一个具体的阈值来代替inf
            lengthsTable(axisIndex, angleIndex) = lengthsTable(axisIndex, angleIndex) + gb(i).segLength;
        end
    end
end
