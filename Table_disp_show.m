function Table_disp_show(gb,axes,angles,lengthsTable,excelname,sheetname)

%% 3.2 将结果转换为table以便更好的展示和使用

% 将结果转换为table以便更好的展示和使用
axisNames = {'[001]','[011]','[111]','[012]','[112]','[122]','[013]','[113]'};
angleNames = arrayfun(@(x) ['Angle', num2str(x)], angles, 'UniformOutput', false);

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

sheetname_lengths = [sheetname,'Lengths'];
sheetname_percent = [sheetname,'Percent'];

% 写入Excel
writetable(lengthsTableUpdated, excelname, 'WriteRowNames', true, 'Sheet', sheetname_lengths);
writetable(percentTableForExcel, excelname, 'WriteRowNames', true, 'Sheet', sheetname_percent);

end