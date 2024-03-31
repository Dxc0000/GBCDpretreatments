clear;clc;

%% 0. 导入并初步处理数据

% 导入数据
filename = "D:\OneDrive - sjtu.edu.cn\课题\表征\EBSD\deformation\data process\500x_crc\AB-front nodeform\xifen-zhengmian-layer0.cpr";
ebsd = EBSD.load(filename,'interface','crc');
EBSD = ebsd;

% 移除不想要的相
removePhaseIds = [0,1,3];
removeInds = ismember(ebsd.phaseId,removePhaseIds);
ebsd(removeInds) = [];

% 计算晶粒
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',5*degree);

% 去除小晶粒
minGrainSize = 2; % 定义最小晶粒大小
ebsd(grains(grains.grainSize < minGrainSize)) = [];

% 重新计算晶粒
[grains, ebsd.grainId] = calcGrains(ebsd,'angle',5*degree);

%% -------------------------------------------------------------------------------

%% 0.1 输出.ang文件，获取信息
% export_ang(ebsd,'D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\AB_front.ang')

%% -------------------------------------------------------------------------------

%% 1. 输出处理后的.ang文件
% % 获取Ni-superalloy的晶体对称性
% cs = ebsd('Ni-superalloy').CS;
% 
% % 复制ebsd对象以保留原始数据（如果需要）
% ebsd_corrected = ebsd;
% 
% % 遍历grains
% for i = 1:length(grains)
%     % 计算当前晶粒的平均取向
%     meanOrientation = grains(i).meanOrientation;
%     
%     % 获取当前晶粒的所有EBSD点
%     currentGrainEbsdIndices = ebsd.grainId == grains(i).id;  % 这是一个逻辑索引
%     
%     % 更新这些EBSD点的取向为平均取向
%     ebsd_corrected(currentGrainEbsdIndices).orientations = repmat(meanOrientation, sum(currentGrainEbsdIndices), 1);
% end
%
% export_ang(ebsd_corrected,'D:\OneDrive\桌面\xifen-zhengmian-layer0.ang')

%% -------------------------------------------------------------------------------

%% 2. 输出处理后的.txt文件，直接用于calc_discrete_dist程序的输出
% smooth这个函数一定要用！如果不用smooth晶界是x，y方向上的折线
grains_smooth = smooth(grains,4);

% 提取晶界
gb = grains_smooth.boundary;

% 初始化矩阵以存储数据
validGbCount = sum(all(gb.ebsdId > 0, 2)); % 只计算两侧都有晶粒的晶界数量
data = zeros(validGbCount, 10); % 初始化存储空间

counter = 0; % 用于跟踪有效晶界段的数量
items = []; % 用于跟踪不符合条件的晶界段序号
for i = 1:length(gb)
    x = gb(i).x; y = gb(i).y;
    if all(gb(i).ebsdId > 0) && numel(unique(x))==numel(x) && numel(unique(y))==numel(y)  % 确保两侧的grainId都大于0
        counter = counter + 1;
        % 提取两侧晶粒的取向（欧拉角）
        ori1 = EBSD(gb(i).ebsdId(1)).orientations.Euler;
        ori2 = EBSD(gb(i).ebsdId(2)).orientations.Euler;
        
        % 提取晶界段的两端坐标
        coords = [gb(i).x(1), gb(i).y(1), gb(i).x(2), gb(i).y(2)];
        
        % 将提取的数据放入data矩阵
        data(counter, :) = [ori1, ori2, coords];

    else
        items = [items i];
    end
end

% 去除不符合条件的gb
gb(items) = [];

% 根据实际收集的数据量调整data大小
data = data(1:counter, :);
roundData = round(data,3);

% 写入数据到.txt，修改文件名为all_segments_Ni.txt
fid = fopen('D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\all_segments_Ni.txt','w');

% 首先，写入指定的文本内容
headerLines = {
    '# File written by extract_gb_traces: version 12/21/2022'
    '# Using data from file: AB_front_000.ang'
    '# data from TSL'
    '# data is on a square grid'
    '# length of step in x direction:    0.650000'
    '# length of step in y direction:    0.650000'
    '# number of rows:  616'
    '# number of columns in odd numbered rows:  769'
    '# number of columns in even numbered rows:  769'
    '# maximal number of pixels in each line segments:    6'
    '# Ni-superalloy Example for gbXstallography'
    '# There are :     543 grains'
    '# There are :    1104 trinodes'
    '# There are :   24604 grain boundary nodes'
    '# There are :    1483 grain pairs with boundaries'
    '# The number of traces subdivided:  114'
    '#     86 of   1483 boundaries do not have 2 trinodes'
    '#      8 of   1483 boundaries have 0 trinode'
    '#      6 of   1483 boundaries have 1 trinode'
    '#   1397 of   1483 boundaries have 2 trinodes'
    '#      0 of   1483 boundaries have 3 trinodes'
    '#     70 of   1483 boundaries have 4 trinodes'
    '#      0 of   1483 boundaries have 5 trinodes'
    '#      1 of   1483 boundaries have 6 trinodes'
    '#      0 of   1483 boundaries have 7 trinodes'
    '#      1 of   1483 boundaries have 8 trinodes'
    '#  phi1    PHI   phi2   phi1    PHI   phi2       x1       y1       x2       y2'
};

for i = 1:numel(headerLines)
    fprintf(fid, '%s\n', headerLines{i});
end

% 然后，遍历矩阵中的每一行，写入数据
for i = 1:size(data, 1)
    % 前6个数据，每个占7个位置
    for j = 1:6
        fprintf(fid, '%7.3f', data(i, j)); % 使用%7.3f确保数字前面有足够的空格，并且数值保留三位小数
    end
    % 后4个数据，每个占9个位置
    for j = 7:10
        fprintf(fid, '%9.3f', data(i, j)); % 使用%9.3f确保数字前面有足够的空格，并且数值保留三位小数
    end
    fprintf(fid, '\n'); % 每行数据结束后换行
end

% 关闭文件
fclose(fid);

%% -------------------------------------------------------------------------------
%% 3. 统计分类晶界段，计算轴角对表格

% 先把轴角对——>旋转矩阵，并投影到基本域中

% 定义晶体对称性
CS = crystalSymmetry('m-3m');

% 定义轴
axes = [Miller(0,0,1,CS), Miller(0,1,1,CS), Miller(-1,1,1,CS), ...
        Miller(0,1,2,CS), Miller(-1,1,2,CS), Miller(-1,2,2,CS), ...
        Miller(0,1,3,CS), Miller(-1,1,3,CS)];

% 定义角度（单位：度）
angles = 5:5:60;

% 初始化轴角对数组
axisAnglePairs = cell(length(axes) * length(angles), 2);

% 生成轴角组合的旋转矩阵
rotMats = cell(length(axes),length(angles));
index = 1;
for i = 1:length(axes)
    for j = 1:length(angles)
        % 创建旋转对象
        rot = rotation('axis', axes(i), 'angle', angles(j)*degree);
        % 转换为旋转矩阵并存储
        rotMats{i,j} = matrix(rot);
        % 存储对应的轴和角度
        axisAnglePairs{index, 1} = axes(i);
        axisAnglePairs{index, 2} = angles(j);
        index = index + 1;
    end
end

% 再把所有晶界段的轴角对——>旋转矩阵，投影到基本域中
mis = project2FundamentalRegion(gb.misorientation);

%% 3.1.1 使用Forbenius方法计算旋转矩阵间的距离

% 初始化Frobenius范数的平方数组
frobeniusSquares = zeros(length(mis), length(rotMats));

% 初始化用于存储最小Frobenius范数的平方的索引数组
minFroIndices = zeros(length(mis), 1);

% 创建一个新的数组来保存对应晶界的最小Frobenius范数的平方的轴角对
alignedAxisAnglePairs = cell(length(mis), 3);

% 初始化一个8x12的表格，用于存储不同轴角对的晶界长度总和
lengthsTable = zeros(length(axes), length(angles));

% 计算每个晶界与参考旋转矩阵的Frobenius范数的平方，并累加长度
for i = 1:length(mis)
    actualRotMat = matrix(mis(i));
    frobeniusSquares = zeros(1, numel(rotMats));
    
    for j = 1:numel(rotMats)
        % 计算Frobenius范数的平方
        frobeniusSquares(j) = norm(actualRotMat - rotMats{j}, 'fro')^2;
    end
    
    [minFrobenius, minIndex] = min(frobeniusSquares);
    
    [axisIndex, angleIndex] = ind2sub([length(axes), length(angles)], minIndex);
    alignedAxisAnglePairs{i, 1} = axes(axisIndex).hkl; % 轴
    alignedAxisAnglePairs{i, 2} = angles(angleIndex); % 角度
    alignedAxisAnglePairs{i, 3} = minFrobenius; % 最小Frobenius范数的平方

    % 如果Frobenius范数的平方小于0.1，累加对应轴角对的晶界长度
    if minFrobenius < 0.1
        lengthsTable(axisIndex, angleIndex) = lengthsTable(axisIndex, angleIndex) + gb(i).segLength;
    end
end


%% 3.1.2 使用旋转角度差方法计算旋转矩阵间的距离

% % 初始化一个8x12的表格，用于存储不同轴角对的晶界长度总和
% lengthsTable = zeros(length(axes), length(angles));
% 
% % 对于每个晶界取向差
% for i = 1:length(mis)
%     actualRotMat = matrix(mis(i)); % 实际取向差的旋转
%     minAngleDiff = inf; % 初始化最小角度差
%     minIndex = 0; % 初始化最小角度差的索引
%     
%     % 对于每个旋转矩阵
%     for j = 1:numel(rotMats)
%         % 计算当前旋转与晶界取向差之间的角度
%         angleDiff = acosd((trace(rotMats{j}' * actualRotMat) - 1) / 2);
%         
%         % 如果角度差小于当前最小值，则更新最小值和索引
%         if angleDiff < minAngleDiff
%             minAngleDiff = angleDiff;
%             minIndex = j;
%         end
%     end
%     
%     % 如果最小角度差小于5度，则累加对应轴角对的晶界长度
%     if minAngleDiff < 5
%         [axisIndex, angleIndex] = ind2sub([length(axes), length(angles)], minIndex);
%         lengthsTable(axisIndex, angleIndex) = lengthsTable(axisIndex, angleIndex) + gb(i).segLength;
%     end
% end

%% 3.2 将结果转换为table以便更好的展示和使用

% 将结果转换为table以便更好的展示和使用
axisNames = arrayfun(@(x) ['Axis_', mat2str(x.hkl)], axes, 'UniformOutput', false);
angleNames = arrayfun(@(x) ['Angle', num2str(x*5)], 1:12, 'UniformOutput', false);

% 在lengthsTable中添加列的和
columnSums = sum(lengthsTable, 1);
lengthsTable = [lengthsTable; columnSums];

% 在lengthsTable中添加行的和
rowSums = sum(lengthsTable, 2);
lengthsTable = [lengthsTable, rowSums];

% 更新轴和角度名称，加入Axis_Sum和Angle_Sum
axisNamesUpdated = [axisNames, {'Axis_Sum'}]; % 假设axisNames已定义
angleNamesUpdated = [angleNames, {'Angle_Sum'}]; % 假设angleNames已定义

% 将更新后的矩阵转换为表格
lengthsTableUpdated = array2table(lengthsTable, 'RowNames', axisNamesUpdated, 'VariableNames', angleNamesUpdated);

% 显示更新后的表格
disp(lengthsTableUpdated);

% 计算总晶界长度
totalLength = sum([gb.segLength]);

% 基于原始长度表（不包括总和行和列）计算百分比
percentTable = lengthsTable / totalLength * 100; % 注意：这里的lengthsTable应是未添加总和前的数值矩阵

% 转换百分比矩阵为表格，准备写入Excel
percentTableForExcel = array2table(percentTable, 'RowNames', axisNamesUpdated, 'VariableNames', angleNamesUpdated);

% 写入Excel
excelname = 'D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\GB_seglength.xlsx';
writetable(lengthsTableUpdated, excelname, 'WriteRowNames', true, 'Sheet', 'Lengths');
writetable(percentTableForExcel, excelname, 'WriteRowNames', true, 'Sheet', 'Percent');

%% -------------------------------------------------------------------------------

%% 4. 绘制取向差分布图
figure('Name','Misorientation')
plotAngleDistribution(mis);
hold on
plotAngleDistribution(CS,CS);
hold off

figurename = 'D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\Mackenzieplot.tif';
print('-dtiff',figurename);
