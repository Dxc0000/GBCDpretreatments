function [DotQuatDiff, minQuatIndices, minQuatAngles, lengthsTable] = DotQuaTable(gb, mis, axes, angles, axisAnglePairs)
    % 计算长度
    Length = length(axes) * length(angles);
    numMis = length(mis);

    % 初始化矩阵和数组
    DotQuatDiff = zeros(numMis, Length);
    minQuatIndices = zeros(numMis, 1);
    minQuatAngles = zeros(numMis, 1);
    lengthsTable = zeros(length(axes), length(angles));
    alignedAxisAnglePairs = cell(numMis, 2);

    % 提取mis中的四元数并转换为矩阵
    quatMis = zeros(numMis, 4);
    for i = 1:numMis
        quatMis(i, :) = [mis(i).a, mis(i).b, mis(i).c, mis(i).d];
    end

    quatAxisAnglePairs = zeros(Length,4);
    for j = 1:Length
        aap = axisAnglePairs{j,3};
        quatAxisAnglePairs(j,:) = [aap.a aap.b aap.c aap.d];
    end

    % 循环计算每个mis和所有axisAnglePairs之间的点积
    for i = 1:numMis
        for j = 1:Length
            % 对当前轴角对四元数与晶界取向差四元数的点积计算夹角
            % 计算点积
            dotProduct = dot(quatMis(i, :), quatAxisAnglePairs(j,:));
            % 计算角度
            DotQuatDiff(i, j) = acosd(dotProduct);
        end

        [minQuatAngles(i), minQuatIndices(i)] = min(DotQuatDiff(i, :));

        [axisIndex, angleIndex] = ind2sub([length(axes), length(angles)], minQuatIndices(i));
        alignedAxisAnglePairs{i, 1} = axes(axisIndex); % 假设axes是一个结构体数组
        alignedAxisAnglePairs{i, 2} = angles(angleIndex); % 假设angles是一个向量

        % 如果四元数之间的夹角小于1度
        if minQuatAngles(i) < 5
            lengthsTable(axisIndex, angleIndex) = lengthsTable(axisIndex, angleIndex) + gb(i).segLength;
        end
    end
end