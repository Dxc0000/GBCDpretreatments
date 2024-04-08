function [RotAngleDiff, minRotIndices, minRotAngles, lengthsTable] = RotAngleDiffTable(gb,mis,rotMats,axes,angles)

    % 初始化Frobenius范数的平方数组
    RotAngleDiff = zeros(length(mis), numel(rotMats));

    % 初始化用于存储最小Frobenius范数的平方的索引数组
    minRotIndices = zeros(length(mis), 1);
    minRotAngles = zeros(length(mis), 1);

    % 初始化一个8x12的表格，用于存储不同轴角对的晶界长度总和
    lengthsTable = zeros(length(axes), length(angles));

    % 创建一个新的数组来保存对应晶界的最小Frobenius范数的平方的轴角对
    alignedAxisAnglePairs = cell(length(mis), 2);

    % 对于每个晶界取向差
    for i = 1:length(mis)
        actualRotMat = matrix(mis(i)); % 实际取向差的旋转
        rotAngleDiff = zeros(1, numel(rotMats));
        
        % 对于每个旋转矩阵
        for j = 1:numel(rotMats)
            % 计算当前旋转与晶界取向差之间的角度
            angleDiff = acosd((trace(rotMats{j}' * actualRotMat) - 1) / 2);
            rotAngleDiff(j) = angleDiff;
        end
        
        [minRotAngle, minRotIndice] = min(rotAngleDiff);
        RotAngleDiff(i,:) = rotAngleDiff;
        minRotIndices(i) = minRotIndice;
        minRotAngles(i) = minRotAngle;

        [axisIndex, angleIndex] = ind2sub([length(axes), length(angles)], minRotIndice);% 线性索引回下标
        alignedAxisAnglePairs{i, 1} = axes(axisIndex).hkl; % 轴
        alignedAxisAnglePairs{i, 2} = angles(angleIndex); % 角度

        % 如果最小角度差小于1度，则累加对应轴角对的晶界长度
        if minRotAngle < 5
            lengthsTable(axisIndex, angleIndex) = lengthsTable(axisIndex, angleIndex) + gb(i).segLength;
        end
    end
end