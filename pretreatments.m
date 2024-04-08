clear;clc;

%% 0. 导入并初步处理数据

% 导入数据
filename = "D:\OneDrive - sjtu.edu.cn\课题\表征\EBSD\deformation\data process\500x_crc\AB-side nodeform\xifen-cemian-300x-layer0.cpr";
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

% smooth这个函数一定要用！如果不用smooth晶界是x，y方向上的折线
grains_smooth = smooth(grains,4);

% 提取晶界
gb = grains_smooth.boundary;

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

%% 2.1 写入数据到.txt，修改文件名为all_segments_Ni.txt
% fid = fopen('D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\all_segments_Ni.txt','w');
% 
% % 首先，写入指定的文本内容
% headerLines = {
%     '# File written by extract_gb_traces: version 12/21/2022'
%     '# Using data from file: AB_front_000.ang'
%     '# data from TSL'
%     '# data is on a square grid'
%     '# length of step in x direction:    0.650000'
%     '# length of step in y direction:    0.650000'
%     '# number of rows:  616'
%     '# number of columns in odd numbered rows:  769'
%     '# number of columns in even numbered rows:  769'
%     '# maximal number of pixels in each line segments:    6'
%     '# Ni-superalloy Example for gbXstallography'
%     '# There are :     543 grains'
%     '# There are :    1104 trinodes'
%     '# There are :   24604 grain boundary nodes'
%     '# There are :    1483 grain pairs with boundaries'
%     '# The number of traces subdivided:  114'
%     '#     86 of   1483 boundaries do not have 2 trinodes'
%     '#      8 of   1483 boundaries have 0 trinode'
%     '#      6 of   1483 boundaries have 1 trinode'
%     '#   1397 of   1483 boundaries have 2 trinodes'
%     '#      0 of   1483 boundaries have 3 trinodes'
%     '#     70 of   1483 boundaries have 4 trinodes'
%     '#      0 of   1483 boundaries have 5 trinodes'
%     '#      1 of   1483 boundaries have 6 trinodes'
%     '#      0 of   1483 boundaries have 7 trinodes'
%     '#      1 of   1483 boundaries have 8 trinodes'
%     '#  phi1    PHI   phi2   phi1    PHI   phi2       x1       y1       x2       y2'
% };
% 
% for i = 1:numel(headerLines)
%     fprintf(fid, '%s\n', headerLines{i});
% end
% 
% % 然后，遍历矩阵中的每一行，写入数据
% for i = 1:size(data, 1)
%     % 前6个数据，每个占7个位置
%     for j = 1:6
%         fprintf(fid, '%7.3f', data(i, j)); % 使用%7.3f确保数字前面有足够的空格，并且数值保留三位小数
%     end
%     % 后4个数据，每个占9个位置
%     for j = 7:10
%         fprintf(fid, '%9.3f', data(i, j)); % 使用%9.3f确保数字前面有足够的空格，并且数值保留三位小数
%     end
%     fprintf(fid, '\n'); % 每行数据结束后换行
% end
% 
% % 关闭文件
% fclose(fid);

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
axisAnglePairs = cell(length(axes) * length(angles), 3);

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
        axisAnglePairs{index, 3} = orientation(rot,CS);
        index = index + 1;
    end
end

% 再把所有晶界段的轴角对——>旋转矩阵，投影到基本域中
mis = project2FundamentalRegion(gb.misorientation);

% 写入Excel
excelname = 'D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\GB_seglength.xlsx';

%% 3.1 使用Forbenius方法计算旋转矩阵间的距离

[frobeniusSquares, minFroIndices, minFroSquares, lengthsTable] = FrobeniusTable(gb,mis,rotMats,axes,angles);

%% 3.1.x 绘制表格
sheetname = 'Frobenius';
Table_disp_show(gb,axes,angles,lengthsTable,excelname,sheetname);

%% 3.2 使用旋转角度差方法计算旋转矩阵间的距离

[RotAngleDiff, minRotIndices, minRotAngles, lengthsTable] = RotAngleDiffTable(gb,mis,rotMats,axes,angles);

%% 3.2.x 绘制表格
sheetname = 'RotAngleDiff';
Table_disp_show(gb,axes,angles,lengthsTable,excelname,sheetname);

%% 3.3 使用最小旋转角度(四元数点积)方法计算旋转矩阵间的距离

[DotQuatDiff, minQuatIndices, minQuatAngles, lengthsTable] = DotQuaTable(gb,mis,axes,angles,axisAnglePairs);

%% 3.3.x 绘制表格
sheetname = 'DotQuatDiff';
Table_disp_show(gb,axes,angles,lengthsTable,excelname,sheetname);

%% 3.4 使用四元数欧式距离法计算旋转矩阵间的距离

[NormQuatDiff, minQuatnIndices, minQuatnAngles, lengthsTable] = NormQuaTable(gb,mis,axes,angles,axisAnglePairs);

%% 3.4.x 绘制表格
sheetname = 'NormQuatDiff';
Table_disp_show(gb,axes,angles,lengthsTable,excelname,sheetname);

%% -------------------------------------------------------------------------------

%% 4. 绘制取向差分布图
figure('Name','Misorientation')
plotAngleDistribution(mis);
hold on
plotAngleDistribution(CS,CS);
hold off

figurename = 'D:\OneDrive\桌面\GSRohrer_programs_transfer_docs\Mackenzieplot.tif';
print('-dtiff',figurename);

close all
