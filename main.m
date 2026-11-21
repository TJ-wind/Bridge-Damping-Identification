%% main.m - 双通道贝叶斯模态辨识、物理参数重构与完整绘图
% main.m 与 main_function.m 放在同一文件夹。
% 在下方填写数据导入语句，运行后先选择两路共同起点，再分别选择终点。
% 图中会比较：原始全程、选中时间段、拟合曲线、残差、频谱、模态参数。

clc;
close all;

%% ===================== 数据输入区（从这里开始） =====================
% 直接在这里输入你要使用的两个通道。

%
% 例 2：如果你已经有变量 beam 和 damper，则改成：
% structureSignal = beam;
% exciterSignal   = damper;
%
% 例 3：如果数据在 MAT 文件里，则改成：
Sdata = load('../data/Bridge2_20260129/test3_190mm_50Hz.mat');
structureSignal = Sdata.beam - mean(Sdata.beam);
exciterSignal   = Sdata.damper - mean(Sdata.damper);
%
fs = 50;  % 采样频率（Hz），修改成真实值

% 物理重构配置（以下对应上面的 Bridge2 加速度数据；更换数据时一并修改）。
physicalOptions = struct();
physicalOptions.m2 = 5000;                 % 已知振子质量，kg
physicalOptions.hasInerter = false;        % Bridge2=false；C240/H240=true
physicalOptions.useKnownAttachment = true; % 有惯容时，固定已知调谐条件
physicalOptions.bareExciterFrequencyHz = 0.5;
physicalOptions.tunedExciterFrequencyHz = 0.3333; % C240=0.2300；H240=0.3333
physicalOptions.responseOrders = [2,2];    % 0位移、1速度、2加速度；Bridge3=[1,1]
physicalOptions.inputScaleToSI = [1,1];    % 原输入转 SI 单位的带符号倍率
% 其余先验和搜索设置见 main_function.m 的“物理重构设置区”。
%% ================================================================

if ~exist('structureSignal','var') || isempty(structureSignal)
    error('请在 main.m 的数据输入区设置 structureSignal，例如 importdata(''beam.txt'')。');
end
if ~exist('exciterSignal','var') || isempty(exciterSignal)
    error('请在 main.m 的数据输入区设置 exciterSignal，例如 importdata(''damper.txt'')。');
end
if exist('main_function','file') ~= 2
    error('找不到 main_function.m。请将它与 main.m 放在同一个文件夹。');
end

fprintf('\n接下来请先选择两路共同起点，再分别选择结构和激振器通道的终点。\n');
result = main_function(structureSignal, exciterSignal, fs, physicalOptions);
physicalFit = result.physicalFit; % 质量、刚度、阻尼、惯容及重构误差，直接留在工作区
if ~isfield(result,'selection')
    error('main_function 缺少选区信息。请下载与此 main.m 配套的新版 main_function.m。');
end
sel = result.selection;

%% 图 1：完整原始信号 + 选区内拟合曲线
figure('Color','w','Name','原始全程与选区拟合');
subplot(2,1,1);
plot(sel.structureTime,sel.structureOriginal,'Color',[0.55 0.55 0.55], ...
    'DisplayName','结构原始全程'); hold on;
plot(sel.structureSelectedTime,result.channels(1).fittedWithOffset,'r-', ...
    'LineWidth',1.4,'DisplayName','选区内拟合');
xline(sel.structureSelectedTime(1),'g--','开始');
xline(sel.structureSelectedTime(end),'m--','结束');
grid on; xlabel('时间（秒）'); ylabel('响应');
title('结构通道：全程原始信号与选区拟合'); legend('Location','best'); hold off;

subplot(2,1,2);
plot(sel.exciterTime,sel.exciterOriginal,'Color',[0.55 0.55 0.55], ...
    'DisplayName','激振器原始全程'); hold on;
plot(sel.exciterSelectedTime,result.channels(2).fittedWithOffset,'r-', ...
    'LineWidth',1.4,'DisplayName','选区内拟合');
xline(sel.exciterSelectedTime(1),'g--','开始');
xline(sel.exciterSelectedTime(end),'m--','结束');
grid on; xlabel('时间（秒）'); ylabel('响应');
title('激振器通道：全程原始信号与选区拟合'); legend('Location','best'); hold off;

%% 图 2：选中时间段内原始信号 vs 拟合（含原始偏置）
figure('Color','w','Name','选区内原始信号与拟合');
subplot(2,1,1);
plot(sel.structureSelectedTime,sel.structureSelected,'b-','DisplayName','选区原始信号'); hold on;
plot(sel.structureSelectedTime,result.channels(1).fittedWithOffset,'r-','LineWidth',1.3, ...
    'DisplayName','拟合曲线（含偏置）');
grid on; xlabel('原始记录时间（秒）'); ylabel('响应');
title(sprintf('结构通道：原始 vs 拟合，时域残差 %.3f%%',100*result.channels(1).relativeTimeError));
legend('Location','best'); hold off;

subplot(2,1,2);
plot(sel.exciterSelectedTime,sel.exciterSelected,'b-','DisplayName','选区原始信号'); hold on;
plot(sel.exciterSelectedTime,result.channels(2).fittedWithOffset,'r-','LineWidth',1.3, ...
    'DisplayName','拟合曲线（含偏置）');
grid on; xlabel('原始记录时间（秒）'); ylabel('响应');
title(sprintf('激振器通道：原始 vs 拟合，时域残差 %.3f%%',100*result.channels(2).relativeTimeError));
legend('Location','best'); hold off;

%% 图 3：去均值后的拟合曲线
figure('Color','w','Name','去均值信号拟合');
for j = 1:2
    c = result.channels(j);
    subplot(2,1,j);
    plot(c.relativeTime,c.observed,'Color',[0.15 0.4 0.75],'DisplayName','观测（去均值）'); hold on;
    plot(c.relativeTime,c.fitted,'r-','LineWidth',1.2,'DisplayName','模型拟合（去均值）');
    grid on; xlabel('选区内相对时间（秒）'); ylabel('响应');
    if j==1, title('结构通道：去均值拟合'); else, title('激振器通道：去均值拟合'); end
    legend('Location','best'); hold off;
end

%% 图 4：时域残差
figure('Color','w','Name','时域残差');
for j = 1:2
    c = result.channels(j);
    residual = c.observed - c.fitted;
    subplot(2,1,j);
    plot(c.relativeTime,residual,'k-'); hold on; yline(0,'r--'); hold off;
    grid on; xlabel('选区内相对时间（秒）'); ylabel('观测 - 拟合');
    title(sprintf('通道 %d 残差，残差 RMS = %.5g',j,sqrt(mean(residual.^2))));
end

%% 图 5：观测频谱 vs 模型频谱
figure('Color','w','Name','频谱拟合对比');
for j = 1:2
    c = result.channels(j);
    subplot(1,2,j);
    semilogy(c.frequencyHz,max(abs(c.spectrum),eps),'.-', ...
        'Color',[0.15 0.4 0.75],'DisplayName','观测频谱'); hold on;
    semilogy(c.frequencyHz,max(abs(c.fittedSpectrum),eps),'-', ...
        'Color',[0.85 0.2 0.15],'LineWidth',1.2,'DisplayName','模型频谱');
    xline(result.frequencyHz(1),'k--','第1阶');
    xline(result.frequencyHz(2),'m--','第2阶');
    grid on; xlabel('频率（Hz）'); ylabel('频谱幅值（对数坐标）');
    if j==1, title('结构通道频谱'); else, title('激振器通道频谱'); end
    legend('Location','best'); hold off;
end

%% 图 6：频谱幅值残差
figure('Color','w','Name','频谱幅值残差');
for j = 1:2
    c = result.channels(j);
    subplot(2,1,j);
    plot(c.frequencyHz,abs(c.spectrum)-abs(c.fittedSpectrum),'b-'); hold on;
    yline(0,'k--'); hold off; grid on;
    xlabel('频率（Hz）'); ylabel('观测幅值 - 拟合幅值');
    if j==1, title('结构通道频谱残差'); else, title('激振器通道频谱残差'); end
end

%% 图 7：残差分布
figure('Color','w','Name','残差分布');
for j = 1:2
    c = result.channels(j);
    residual = c.observed-c.fitted;
    subplot(1,2,j);
    histogram(residual,30);
    grid on; xlabel('残差'); ylabel('点数');
    if j==1, title('结构通道残差分布'); else, title('激振器通道残差分布'); end
end

%% 图 8：模态频率、阻尼比和振型
figure('Color','w','Name','模态参数与振型');
subplot(2,2,1);
bar(1:2,result.frequencyHz);
set(gca,'XTick',1:2,'XTickLabel',{'第1阶','第2阶'});
ylabel('固有频率（Hz）'); title('固有频率'); grid on;

subplot(2,2,2);
bar(1:2,result.dampingPercent);
set(gca,'XTick',1:2,'XTickLabel',{'第1阶','第2阶'});
ylabel('阻尼比（%）'); title('阻尼百分比'); grid on;

subplot(2,2,3);
bar(abs(result.modeShapes));
set(gca,'XTick',1:2,'XTickLabel',{'结构','激振器'});
ylabel('归一化振型幅值'); title('复数模态振型幅值'); legend('第1阶','第2阶'); grid on;

subplot(2,2,4);
bar(angle(result.modeShapes)*180/pi);
set(gca,'XTick',1:2,'XTickLabel',{'结构','激振器'});
ylabel('相位（度）'); title('复数模态振型相位'); legend('第1阶','第2阶'); grid on;

%% 保存 main.m 额外生成的所有图和摘要
outRoot = fullfile(fileparts(mfilename('fullpath')),'main_plots');
if exist(outRoot,'dir') ~= 7, mkdir(outRoot); end
outFolder = fullfile(outRoot,['run_' datestr(now,'yyyymmdd_HHMMSS')]);
mkdir(outFolder);
figs = findall(0,'Type','figure');
for k = 1:numel(figs)
    try
        saveas(figs(k),fullfile(outFolder,sprintf('figure_%02d.png',k)));
    catch
        % 单张图保存失败时继续保存其余图
    end
end
save(fullfile(outFolder,'result_with_selection.mat'),'result','physicalFit','sel','fs');
writetable(result.modalTable,fullfile(outFolder,'modal_parameters.csv'));
writetable(physicalFit.physicalParameterTable,fullfile(outFolder,'physical_parameters.csv'));
writetable(physicalFit.parameterTable,fullfile(outFolder,'physical_search_parameters.csv'));
writetable(physicalFit.qualityTable,fullfile(outFolder,'physical_fit_quality.csv'));
fprintf('\nmain.m 的附加图和结果已保存到：\n%s\n',outFolder);
fprintf('已自动导出物理参数、搜索参数和物理重构残差三张 CSV 表。\n');
fprintf('注意：模型假设两个稳定欠阻尼模态；图形拟合好并不自动证明模型完全适合实测数据。\n');
