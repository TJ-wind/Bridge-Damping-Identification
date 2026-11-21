function result = main_function(structureSignal, exciterSignal, fs, physicalOptions)
%MAIN_FUNCTION 双通道贝叶斯模态辨识与物理参数重构入口
%
% 用法A：信号已经在工作区中：
%   result = main_function(beam, damper, 50);
%   result = main_function(beam, damper, 50, physicalOptions);
% physicalOptions 可覆盖下方“物理重构设置区”的同名字段；默认自动重构。
% 例如 Bridge2 加速度记录：struct('hasInerter',false,'responseOrders',[2,2])。
% 仅进行原谱辨识时，可传 struct('enabled',false)。
%
% 用法B：在下面“数据输入区”取消注释并改成自己的文件/变量，
%         然后运行 main_function：
%   structureSignal = importdata('beam.txt');
%   exciterSignal   = importdata('damper.txt');
% 或：
%   S = load('mydata.mat');
%   structureSignal = S.beam;
%   exciterSignal   = S.damper;
%
% 本文件包含数据输入、鼠标选区、原始辨识算法、物理重构入口和结果保存。
% 不需要额外的 bayesian_spectrum_identification.m 或绘图辅助文件。
% 物理重构复用同目录 refit_record5_unequal_free_k2.m 和 +physical_refit_support，
% 需要 Optimization Toolbox；只计算并导出表格，不新增物理重构图窗。
% 验证脚本和噪声试验未合并，因为它们不是日常实测辨识的必要步骤。
%
% ====================== 数据输入区 ======================
% 如果使用 importdata/load，请把对应语句放在下面 nargin==0 的代码位置。
% ========================================================

    narginchk(0,4);
    if nargin < 4, physicalOptions = struct(); end
    physicalOptions = local_physical_options(physicalOptions);

    if nargin == 0
        % 先初始化为空，再在下面取消注释并修改你需要的输入方式。
        structureSignal = [];
        exciterSignal = [];
        fs = 50; % 修改成实际采样频率（Hz）

        % 方式 A：从文本文件导入（取消注释并修改文件名）
        % structureSignal = importdata('beam.txt');
        % exciterSignal   = importdata('damper.txt');

        % 方式 B：从 MAT 文件导入（取消注释并修改变量名/文件名）
        % S = load('mydata.mat');
        % structureSignal = S.beam;
        % exciterSignal   = S.damper;
    elseif nargin == 1 || nargin == 2
        error('请同时输入两个信号和采样频率，例如 main_function(beam, damper, 50)。');
    end

    structureSignal = unwrap_importdata(structureSignal, '结构通道');
    exciterSignal = unwrap_importdata(exciterSignal, '激振器/振子通道');
    if isempty(structureSignal) || isempty(exciterSignal)
        error(['没有输入信号。请填写函数中的数据导入语句，' ...
            '或者调用 main_function(beam, damper, 50)。']);
    end
    if ~isnumeric(structureSignal) || ~isvector(structureSignal) || ~isreal(structureSignal) || any(~isfinite(structureSignal(:)))
        error('结构通道必须是有限的实数向量。');
    end
    if ~isnumeric(exciterSignal) || ~isvector(exciterSignal) || ~isreal(exciterSignal) || any(~isfinite(exciterSignal(:)))
        error('激振器通道必须是有限的实数向量。');
    end
    validateattributes(fs, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'fs');
    structureSignal = double(structureSignal(:));
    exciterSignal = double(exciterSignal(:));

    fprintf('\n===== 双通道贝叶斯模态辨识 =====\n');
    fprintf('结构通道：%d 点，约 %.4f 秒\n', numel(structureSignal), numel(structureSignal)/fs);
    fprintf('激振器通道：%d 点，约 %.4f 秒\n', numel(exciterSignal), numel(exciterSignal)/fs);
    fprintf('采样频率：%.8g Hz\n', fs);
    fprintf('请先在结构通道选择两路共同起点，再分别选择两路的终点。\n');
    fprintf('两路原始信号的第一个采样点须对应同一实际时刻。\n');

    tB = (0:numel(structureSignal)-1)'/fs;
    tE = (0:numel(exciterSignal)-1)'/fs;
    commonLength = min(numel(structureSignal),numel(exciterSignal));
    if commonLength < 32
        error('两路记录的重叠区间至少需要 32 个采样点。');
    end
    latestStartTime = (commonLength-32)/fs;

    pickFig = figure('Color','w','Name','选择自由衰减区间');
    axB = subplot(2,1,1);
    plot(axB,tB,structureSignal,'b-'); grid(axB,'on');
    xlabel(axB,'时间（秒）'); ylabel(axB,'响应');
    title(axB,{'第 1 步：选择两路共同起点','请在结构通道曲线上点击一个位置'});
    axE = subplot(2,1,2);
    plot(axE,tE,exciterSignal,'Color',[0.85 0.33 0.1]); grid(axE,'on');
    xlabel(axE,'时间（秒）'); ylabel(axE,'响应');
    title(axE,{'激振器/振子通道','与上方通道共用同一个实际起点'});
    drawnow;

    while true
        axes(axB); %#ok<LAXES>
        [xStart,~] = ginput(1);
        if numel(xStart) ~= 1
            if isvalid(pickFig), close(pickFig); end
            error('未完成共同起点选择，请重新运行。');
        end
        if ~isfinite(xStart) || xStart < 0 || xStart > tB(commonLength)
            title(axB,{'共同起点无效，请重新点击', ...
                sprintf('共同起点须在 0 至 %.4f 秒之间，且两路各保留至少 32 点',latestStartTime)});
            warning('共同起点必须位于两路重叠区间内，并为两路各保留至少 32 个采样点。');
            continue;
        end
        commonStartIndex = round(xStart*fs)+1;
        if commonStartIndex < 1 || commonStartIndex+31 > commonLength
            title(axB,{'共同起点无效，请重新点击','两路选区各至少需要 32 个采样点'});
            warning('共同起点对应的采样点无效，请重新选择。');
            continue;
        end
        break;
    end
    commonStartTime = tB(commonStartIndex);
    hold(axB,'on'); xline(axB,commonStartTime,'g--','共同起点'); hold(axB,'off');
    hold(axE,'on'); xline(axE,commonStartTime,'g--','共同起点'); hold(axE,'off');
    fprintf('共同起点：%.6g 秒\n',commonStartTime);

    title(axB,{'第 2 步：选择结构通道终点','请在共同起点之后点击结构通道曲线'});
    title(axE,{'激振器/振子通道','请先选择上方结构通道终点'});
    drawnow;
    while true
        axes(axB); %#ok<LAXES>
        [xBEnd,~] = ginput(1);
        if numel(xBEnd) ~= 1
            if isvalid(pickFig), close(pickFig); end
            error('未完成结构通道终点选择，请重新运行。');
        end
        if ~isfinite(xBEnd) || xBEnd < commonStartTime || xBEnd > tB(end)
            title(axB,{'结构通道终点无效，请重新点击','终点须位于共同起点之后及记录范围内'});
            warning('结构通道终点须位于共同起点之后及记录范围内。');
            continue;
        end
        structureEndIndex = round(xBEnd*fs)+1;
        if structureEndIndex-commonStartIndex+1 < 32
            title(axB,{'结构通道选区过短，请重新选择终点','选区至少需要 32 个采样点'});
            warning('结构通道选区至少需要 32 个采样点。');
            continue;
        end
        break;
    end
    idxB = [commonStartIndex;structureEndIndex];
    hold(axB,'on'); xline(axB,tB(idxB(2)),'r--','结束'); hold(axB,'off');
    title(axB,'结构通道选区已确定');

    title(axE,{'第 3 步：选择激振器/振子通道终点','请在共同起点之后点击激振器/振子曲线'});
    drawnow;
    while true
        axes(axE); %#ok<LAXES>
        [xEEnd,~] = ginput(1);
        if numel(xEEnd) ~= 1
            if isvalid(pickFig), close(pickFig); end
            error('未完成激振器通道终点选择，请重新运行。');
        end
        if ~isfinite(xEEnd) || xEEnd < commonStartTime || xEEnd > tE(end)
            title(axE,{'激振器通道终点无效，请重新点击','终点须位于共同起点之后及记录范围内'});
            warning('激振器通道终点须位于共同起点之后及记录范围内。');
            continue;
        end
        exciterEndIndex = round(xEEnd*fs)+1;
        if exciterEndIndex-commonStartIndex+1 < 32
            title(axE,{'激振器通道选区过短，请重新选择终点','选区至少需要 32 个采样点'});
            warning('激振器通道选区至少需要 32 个采样点。');
            continue;
        end
        break;
    end
    idxE = [commonStartIndex;exciterEndIndex];
    hold(axE,'on'); xline(axE,tE(idxE(2)),'r--','结束'); hold(axE,'off');

    fprintf('结构通道选区：%.6g 到 %.6g 秒\n',tB(idxB(1)),tB(idxB(2)));
    fprintf('激振器通道选区：%.6g 到 %.6g 秒\n',tE(idxE(1)),tE(idxE(2)));
    pause(0.5);
    if isvalid(pickFig), close(pickFig); end

    structureSelected = structureSignal(idxB(1):idxB(2));
    exciterSelected = exciterSignal(idxE(1):idxE(2));
    startTime = 0;
    structureStopTime = numel(structureSelected)/fs;
    exciterStopTime = numel(exciterSelected)/fs;

    fprintf('\n正在辨识，请稍候……\n');
    result = local_bayesian_spectrum_identification( ...
        structureSelected, exciterSelected, startTime, structureStopTime, exciterStopTime, fs);

    % 保存原始信号与手动选区，供 main.m 绘制全程/选区拟合对比图。
    result.selection = struct( ...
        'structureOriginal', structureSignal, 'exciterOriginal', exciterSignal, ...
        'structureTime', tB, 'exciterTime', tE, ...
        'structureIndices', idxB, 'exciterIndices', idxE, ...
        'structureSelected', structureSelected, 'exciterSelected', exciterSelected, ...
        'structureSelectedTime', tB(idxB(1):idxB(2)), ...
        'exciterSelectedTime', tE(idxE(1):idxE(2)));

    fprintf('\n===== 辨识结果 =====\n');
    fprintf('第 1 阶固有频率：%.8g Hz\n', result.frequencyHz(1));
    fprintf('第 2 阶固有频率：%.8g Hz\n', result.frequencyHz(2));
    fprintf('第 1 阶阻尼比：%.8g（%.6g%%）\n', result.dampingRatio(1), result.dampingPercent(1));
    fprintf('第 2 阶阻尼比：%.8g（%.6g%%）\n', result.dampingRatio(2), result.dampingPercent(2));
    fprintf('结构通道时域残差：%.4g%%\n', 100*result.channels(1).relativeTimeError);
    fprintf('激振器通道时域残差：%.4g%%\n', 100*result.channels(2).relativeTimeError);
    fprintf('结构通道频谱残差：%.4g%%\n', 100*result.channels(1).relativeSpectrumError);
    fprintf('激振器通道频谱残差：%.4g%%\n', 100*result.channels(2).relativeSpectrumError);
    fprintf('不确定性诊断：%s\n', result.diagnostics.status);
    if result.diagnostics.intervalsAccepted
        disp('95% 局部近似区间（行：f1、f2、阻尼比1、阻尼比2；列：下限、上限）：');
        disp(result.diagnostics.ci95);
    else
        fprintf('局部区间检查未通过；不要将区间当作有效置信区间使用。\n');
    end
    disp('模态参数表：'); disp(result.modalTable);
    disp('复数模态振型（行：结构/激振器；列：模态）：'); disp(result.modeShapes);

    if physicalOptions.enabled
        fprintf('\n正在重构结构模态质量与物理参数，请稍候……\n');
        [result.physicalFit,physicalOptions] = local_physical_refit(result,physicalOptions);
        fprintf('\n===== 结构模态质量与物理参数重构 =====\n');
        fprintf('结构模态质量 mB = %.10g kg\n',result.physicalFit.model.m1);
        fprintf('物理反演诊断：%s\n',result.physicalFit.status);
        disp('物理参数表：'); disp(result.physicalFit.physicalParameterTable);
        disp('物理搜索参数表（含结构/振子频率、阻尼比及搜索范围）：');
        disp(result.physicalFit.parameterTable);
        disp('物理重构残差表（Percent 列为百分数）：');
        disp(result.physicalFit.qualityTable);
    end
    result.physicalOptions = physicalOptions;

    fig = figure('Color','w','Name','贝叶斯双模态辨识结果');
    for j = 1:2
        c = result.channels(j);
        subplot(2,2,j);
        plot(c.time,c.observed,'Color',[0.6 0.6 0.6]); hold on;
        plot(c.time,c.fitted,'b','LineWidth',1.2); hold off; grid on;
        xlabel('时间（秒）'); ylabel('去均值响应');
        if j==1, title('结构通道时域拟合'); else, title('激振器通道时域拟合'); end
        legend('输入信号','两阶模型拟合','Location','best');

        subplot(2,2,j+2);
        plot(c.frequencyHz,abs(c.spectrum),'Color',[0.6 0.6 0.6]); hold on;
        plot(c.frequencyHz,abs(c.fittedSpectrum),'b','LineWidth',1.2);
        xline(result.frequencyHz(1),'r--','第1阶');
        xline(result.frequencyHz(2),'m--','第2阶'); hold off; grid on;
        xlabel('频率（Hz）'); ylabel('频谱幅值');
        xlim([0,min(fs/2,max(1,2*max(result.frequencyHz)))]);
        if j==1, title('结构通道频谱拟合'); else, title('激振器通道频谱拟合'); end
        legend('输入频谱','拟合频谱','Location','best');
    end

    outRoot = fullfile(fileparts(mfilename('fullpath')),'main_function_results');
    if exist(outRoot,'dir') ~= 7, mkdir(outRoot); end
    outFolder = fullfile(outRoot,['result_' datestr(now,'yyyymmdd_HHMMSS')]);
    mkdir(outFolder);
    result.outputFolder = outFolder;
    save(fullfile(outFolder,'result.mat'),'result','fs','startTime','structureStopTime','exciterStopTime');
    writetable(result.modalTable,fullfile(outFolder,'modal_parameters.csv'));
    if physicalOptions.enabled
        physicalFit = result.physicalFit;
        save(fullfile(outFolder,'physical_fit.mat'),'physicalFit','physicalOptions');
        writetable(physicalFit.physicalParameterTable,fullfile(outFolder,'physical_parameters.csv'));
        writetable(physicalFit.parameterTable,fullfile(outFolder,'physical_search_parameters.csv'));
        writetable(physicalFit.qualityTable,fullfile(outFolder,'physical_fit_quality.csv'));
    end
    saveas(fig,fullfile(outFolder,'fit_plot.png'));
    fprintf('\n结果已保存到：\n%s\n',outFolder);
    fprintf('文件包括 result.mat、modal_parameters.csv、fit_plot.png。\n');
    if physicalOptions.enabled
        fprintf('物理重构包括 physical_fit.mat、physical_parameters.csv、physical_search_parameters.csv、physical_fit_quality.csv。\n');
    end
end

function options = local_physical_options(options)
% ====================== 物理重构设置区 ======================
% 默认沿用 main_physical_refit.m；main.m 可通过第4参数传入当前工况。
% Bridge2：hasInerter=false、responseOrders=[2,2]。
% Bridge3：responseOrders=[1,1]；C240 调谐频率0.2300，H240为0.3333 Hz。
    defaults = struct( ...
        'enabled',true, ...
        'm2',5000, ...                         % 已知振子质量，kg
        'hasInerter',true, ...
        'useKnownAttachment',true, ...
        'bareExciterFrequencyHz',0.5, ...      % 裸振子固定底座频率，Hz
        'tunedExciterFrequencyHz',0.3333, ...  % 设备调谐值，Hz
        'responseOrders',[1,1], ...            % 0位移、1速度、2加速度
        'inputScaleToSI',[1,1], ...            % 带符号的 SI 转换倍率
        'frequencyBandHz',[], ...             % []：由本次谱结果计算
        'priorLower',[], ...
        'priorUpper',[], ...
        'initialParameters',[], ...
        'tau',1e3, ...
        'a0',2, ...
        'b0Normalized',1e-4, ...
        'b0',[], ...                          % []：归一化尺度转为 SI 响应平方
        'stateFrequencyReferenceHz',[], ...
        'statePriorScale',[], ...
        'initialState',[], ...
        'searchSeeds',[202610111,202610112,202610113], ...
        'populationSize',96, ...
        'maxIterations',350, ...
        'maxStallIterations',120, ...
        'localStarts',12, ...
        'localMaxIterations',800);
% ===========================================================
    validateattributes(options,{'struct'},{'scalar'},mfilename,'physicalOptions');
    names = fieldnames(options);
    for k = 1:numel(names)
        if ~isfield(defaults,names{k})
            error('mainFunction:PhysicalOption','未知的物理重构设置：%s。',names{k});
        end
        defaults.(names{k}) = options.(names{k});
    end
    options = defaults;
    for name = {'enabled','hasInerter','useKnownAttachment'}
        value = options.(name{1});
        validateattributes(value,{'logical','numeric'},{'real','scalar','finite','binary'}, ...
            mfilename,['physicalOptions.' name{1}]);
        options.(name{1}) = logical(value);
    end
    if options.enabled
        if exist('refit_record5_unequal_free_k2','file') ~= 2 || ...
                isempty(which('physical_refit_support.objective'))
            error('mainFunction:PhysicalDependency', ...
                '物理重构需要同目录的 refit_record5_unequal_free_k2.m 和 +physical_refit_support 文件夹。');
        end
        if ~license('test','Optimization_Toolbox') || exist('fmincon','file') ~= 2
            error('mainFunction:PhysicalToolbox', ...
                '物理重构需要 Optimization Toolbox (fmincon)。仅做谱辨识可设置 physicalOptions.enabled=false。');
        end
    end
end

function [physicalFit,options] = local_physical_refit(result,options)
% 直接使用本次选段及真实起点，不再导入数据，不调用脚本的绘图部分。
    bridgeResponse = result.selection.structureSelected(:);
    exciterResponse = result.selection.exciterSelected(:);
    fs = result.samplingFrequency;
    startTimes = [result.selection.structureSelectedTime(1), ...
                  result.selection.exciterSelectedTime(1)];
    fmin = min(result.frequencyHz); fmax = max(result.frequencyHz);
    if isempty(options.frequencyBandHz)
        options.frequencyBandHz = [0.3*fmin,min(1.8*fmax,0.49*fs)];
    end
    if isempty(options.priorLower)
        options.priorLower = [1e-3,0.5*fmin,1e-5,1e-5,0,0.01,0.5*fmin];
    end
    if isempty(options.priorUpper)
        options.priorUpper = [1e8,2*fmax,0.3,0.5,0.9999,1000,1.5*fmax];
    end
    if isempty(options.initialParameters)
        options.initialParameters = [30,fmax,0.01,0.05,0.5,1,fmin; ...
                                     30,fmin,0.01,0.05,0.5,1,fmax];
    end
    validateattributes(options.priorLower,{'numeric'},{'real','vector','numel',7,'finite'});
    validateattributes(options.priorUpper,{'numeric'},{'real','vector','numel',7,'finite'});
    validateattributes(options.initialParameters,{'numeric'},{'real','2d','ncols',7,'finite','nonempty'});
    options.priorLower = options.priorLower(:).';
    options.priorUpper = options.priorUpper(:).';
    if ~options.hasInerter
        options.priorLower(5)=0; options.priorUpper(5)=0;
        options.initialParameters(:,5)=0;
    elseif options.useKnownAttachment
        validateattributes(options.bareExciterFrequencyHz,{'numeric'},{'real','scalar','finite','positive'});
        validateattributes(options.tunedExciterFrequencyHz,{'numeric'},{'real','scalar','finite','positive'});
        assert(options.tunedExciterFrequencyHz<=options.bareExciterFrequencyHz, ...
            'mainFunction:Tuning','已知调谐频率不能高于裸振子频率。');
        knownBeta = 1-(options.tunedExciterFrequencyHz/options.bareExciterFrequencyHz)^2;
        options.priorLower([2,5])=[options.tunedExciterFrequencyHz,knownBeta];
        options.priorUpper([2,5])=options.priorLower([2,5]);
        options.initialParameters(:,2)=options.tunedExciterFrequencyHz;
        options.initialParameters(:,5)=knownBeta;
    end
    validateattributes(options.responseOrders,{'numeric'},{'real','vector','numel',2,'integer','>=',0,'<=',2});
    validateattributes(options.inputScaleToSI,{'numeric'},{'real','vector','numel',2,'finite','nonzero'});
    if isempty(options.b0)
        validateattributes(options.b0Normalized,{'numeric'},{'real','scalar','finite','positive'});
        bridgeSI = bridgeResponse*options.inputScaleToSI(1);
        options.b0 = options.b0Normalized*mean((bridgeSI-mean(bridgeSI)).^2);
    end
    if isempty(options.stateFrequencyReferenceHz)
        options.stateFrequencyReferenceHz = fmin;
    end
    if isempty(options.statePriorScale)
        validateattributes(options.stateFrequencyReferenceHz,{'numeric'},{'real','scalar','finite','positive'});
        w = 2*pi*options.stateFrequencyReferenceHz;
        p = options.responseOrders(1);
        options.statePriorScale = [w^p,w^p,w^(p-1),w^(p-1)];
    end
    physicalFit = refit_record5_unequal_free_k2( ...
        bridgeResponse,exciterResponse,fs,startTimes,options.m2,options.hasInerter, ...
        options.responseOrders,options.inputScaleToSI,options.frequencyBandHz, ...
        options.priorLower,options.priorUpper,options.tau,options.a0,options.b0, ...
        options.statePriorScale,options.initialParameters,options.initialState, ...
        options.searchSeeds,options.populationSize,options.maxIterations, ...
        options.maxStallIterations,options.localStarts,options.localMaxIterations);
    timeError = physicalFit.centeredTimeNRMSE(:);
    spectrumError = physicalFit.selectedNRMSE(:);
    physicalFit.qualityTable = table(["structure";"exciter"], ...
        timeError,100*timeError,spectrumError,100*spectrumError, ...
        'VariableNames',{'Channel','TimeNRMSE','TimeNRMSEPercent', ...
        'SelectedSpectrumNRMSE','SelectedSpectrumNRMSEPercent'});
end

function value = unwrap_importdata(value, channelName)
% 处理 importdata 返回的数值数组或含 data 字段的结构体。
    if isstruct(value)
        if isfield(value,'data') && isnumeric(value.data)
            value = value.data;
        else
            error('%s 的 importdata 返回了结构体，但没有可用的数值 data 字段。请手动提取数据。',channelName);
        end
    end
    if isnumeric(value) && ismatrix(value) && size(value,2)>1
        error('%s 是多列矩阵，请明确指定其中代表信号的那一列。',channelName);
    end
end


% ===== 原始贝叶斯辨识核心算法，已整合为本文件的局部函数 =====
function result = local_bayesian_spectrum_identification(structuralSignal,exciterSignal, ...
    startTime,structuralStopTime,exciterStopTime,fs)
%BAYESIAN_SPECTRUM_IDENTIFICATION Two-mode Bayesian free-decay identification.
% R = bayesian_spectrum_identification(yB,yE,t0,tB,tE,fs)
% has exactly six inputs and needs only base MATLAB. This file is standalone.
%
% INPUTS
%   yB,yE   Real vectors of synchronized structural and exciter MOTION, both
%           measured in the same response type (displacement, velocity, or
%           acceleration). Sample 1 of each vector is at time zero. Vectors
%           may have different lengths. The exciter channel is not force.
%   t0      Common start of force-free decay, in seconds from sample 1.
%   tB,tE   Separate exclusive stopping times, in seconds. Samples satisfying
%           t0 <= (index-1)/fs < stop are retained. Off-grid limits therefore
%           select the next sample; limits within roundoff of a grid point
%           are snapped to that point. Either channel may end first.
%   fs      Common sampling frequency in Hz.
%
% OUTPUTS (modes are ordered by increasing coupled natural frequency)
%   frequencyHz, dampedFrequencyHz, dampingRatio, dampingPercent, poles
%   residues       2-by-2 complex response residues: modes-by-channels.
%                  Referenced to the FIRST RETAINED sample in both channels.
%   modeShapes     2-by-2 complex shapes: channels-by-modes, normalized by
%                  the largest component; shapeReferenceChannel gives it.
%   channels(j)    Observed/fitted centered time histories, raw modal mean,
%                  fittedWithOffset, native spectra, indices and fit errors.
%   diagnostics    Local Laplace checks. ci95 is NaN when checks fail;
%                  accepted intervals are local approximations, not MCMC.
%   prior, search  All automatic prior choices and multistart results.
%
% MODEL AND DEFAULTS
%   Exactly TWO stable underdamped conjugate pole pairs are assumed, shared
%   by the channels, with unrestricted complex response residues. Thus
%   non-classical damping is allowed. These are coupled modal parameters,
%   not uncoupled bridge/exciter frequencies or mechanical masses/stiffnesses.
%   Channel calibration and directions must be consistent for mode shapes.
%   No differentiation, integration, tapering or zero padding is performed.
%   Every strictly positive non-Nyquist native DFT bin enters the likelihood.
%   DC is excluded (an unknown constant offset); even-record Nyquist is also
%   excluded to keep all Fourier coordinates identically normalized.
%
%   For each channel: Y = D*w + e; e ~ N(0,v*I),
%   w|v ~ N(0,v*tau^2*I4), v ~ IG(a0,b0), independently across channels.
%   Residues and variance are integrated analytically, including the log
%   determinant term. Poles maximize this marginal posterior. Uniform
%   physical-coordinate priors use 1e-5 < damping < 0.3 and natural frequency
%   max(fs/(100*Nmax),1e-8*fs) < f < (1-1e-8)*fs/2, with f1 < f2.
%   tau=1000, a0=2, b0=1e-4 in channel-RMS-normalized units. In original
%   units b0_j=1e-4*RMS_j^2: this is an explicitly empirical-Bayes default.
%   Separate channel noise variances extend the common-variance modal model
%   in bayesian_spectrum_derivation.tex. These defaults do not reproduce the
%   mechanically constrained field fits or their selected-band priors.
%
% EXAMPLE
%   R = bayesian_spectrum_identification(yB,yE,930,1450,944,50);
%   disp(R.modalTable)
%   plot(R.channels(1).time,R.channels(1).observed); hold on
%   plot(R.channels(1).time,R.channels(1).fitted); hold off
%
% The native DFT convention follows MATLAB fft; optimization uses base
% MATLAB fminsearch (multistart local search, no global-optimum guarantee):
% https://www.mathworks.com/help/matlab/ref/fft.html
% https://www.mathworks.com/help/matlab/ref/fminsearch.html

narginchk(6,6);
clockStart = tic;
validateattributes(structuralSignal,{'numeric'},{'real','vector','finite','nonempty'},mfilename,'structuralSignal');
validateattributes(exciterSignal,{'numeric'},{'real','vector','finite','nonempty'},mfilename,'exciterSignal');
validateattributes(fs,{'numeric'},{'real','scalar','finite','positive'},mfilename,'fs');
validateattributes(startTime,{'numeric'},{'real','scalar','finite','nonnegative'},mfilename,'startTime');
validateattributes(structuralStopTime,{'numeric'},{'real','scalar','finite','positive'},mfilename,'structuralStopTime');
validateattributes(exciterStopTime,{'numeric'},{'real','scalar','finite','positive'},mfilename,'exciterStopTime');
fs=double(fs); startTime=double(startTime);
stops=double([structuralStopTime,exciterStopTime]);
if any(stops<=startTime)
    error('bayesianSpectrum:WindowOrder','Each stopping time must exceed the common starting time.');
end
signals={double(structuralSignal(:)),double(exciterSignal(:))};
channelNames={'structure','exciter'};
data=cell(1,2); channels=struct([]);
first=grid_ceiling(startTime*fs)+1;
for j=1:2
    stopSample=stops(j)*fs;
    if ~isfinite(stopSample) || stopSample>numel(signals{j})+16*eps(max(1,stopSample))
        error('bayesianSpectrum:WindowOutsideSignal', ...
            '%s stop exceeds its available duration, numel(signal)/fs.',channelNames{j});
    end
    last=grid_ceiling(stopSample);
    if first<1 || last>numel(signals{j}) || last-first+1<32
        error('bayesianSpectrum:ShortWindow','Each window must contain at least 32 available samples.');
    end
    indices=(first:last).'; y=signals{j}(indices); N=numel(y);
    removedMean=mean(y); y=y-removedMean;
    scale=norm(y)/sqrt(N);
    if ~isfinite(scale) || scale<=32*eps(max(abs(signals{j}(indices))))
        error('bayesianSpectrum:ConstantSignal','Each retained channel must have nonzero variation.');
    end
    k=(1:floor((N-1)/2)).';
    Z=fft(y/scale)/sqrt(N); Z=Z(k+1);
    if norm(Z)<=32*eps(sqrt(N))
        error('bayesianSpectrum:EmptySpectrum', ...
            'Each channel must contain variation away from DC and Nyquist.');
    end
    data{j}=struct('N',N,'omega',2*pi*k/N,'frequencyHz',k*fs/N, ...
        'Y',sqrt(2)*[real(Z);imag(Z)],'Z',Z,'scale',scale);
    channels(j).name=channelNames{j};
    channels(j).sampleIndices=indices;
    channels(j).sampleCount=N;
    channels(j).time=(indices-1)/fs;
    channels(j).relativeTime=(0:N-1).'/fs;
    channels(j).requestedWindow=[startTime,stops(j)];
    channels(j).effectiveWindow=[(first-1)/fs,last/fs];
    channels(j).observed=y;
    channels(j).removedMean=removedMean;
    channels(j).scale=scale;
    channels(j).frequencyHz=data{j}.frequencyHz;
    channels(j).spectrum=Z*scale;
end
Nmax=max([data{1}.N,data{2}.N]);
lo=[repmat(max(fs/(100*Nmax),fs*1e-8),1,2),1e-5,1e-5];
hi=[repmat((1-1e-8)*fs/2,1,2),.3,.3];
prior=struct('lower',lo,'upper',hi,'tau',1000,'a0',2,'b0Normalized',1e-4, ...
    'b0OriginalUnits',1e-4*[data{1}.scale,data{2}.scale].^2, ...
    'measure','Uniform ordered natural frequencies and damping ratios in physical coordinates', ...
    'scaleChoice','Empirical Bayes: b0 is 1e-4 times each retained centered mean square', ...
    'noiseModel','Independent channel white Gaussian noise and inverse-gamma variance priors');
starts=automatic_starts(data,lo,hi,fs);
objective=@(theta)marginal_objective(theta,data,prior,fs);
options=optimset('Display','off','MaxIter',1600,'MaxFunEvals',2600, ...
    'TolX',1e-8,'TolFun',1e-9);
count=size(starts,1); values=inf(count,1); candidates=nan(count,4);
flags=zeros(count,1); iterations=zeros(count,1);
for s=1:count
    initial=log((starts(s,:)-lo)./(hi-starts(s,:)));
    [u,values(s),flags(s),info]=fminsearch(@(u)objective(decode(u,lo,hi)),initial,options);
    candidates(s,:)=decode(u,lo,hi); iterations(s)=info.iterations;
end
[best,bestIndex]=min(values);
if ~isfinite(best) || best>=1e90
    error('bayesianSpectrum:OptimizationFailed','No finite posterior candidate was found.');
end
theta=candidates(bestIndex,:);
% A final restart with a smaller simplex improves stationarity at the winner.
span=hi-lo; localScale=max(abs(theta),[fs/Nmax,fs/Nmax,.01,.01]);
refineOptions=optimset(options,'MaxFunEvals',4000,'MaxIter',2200,'TolX',1e-9);
[v,refined,flag,info]=fminsearch(@(v)objective(theta+localScale.*v),zeros(1,4),refineOptions);
if refined<=best && isfinite(refined)
    theta=theta+localScale.*v; best=refined;
    candidates(end+1,:)=theta; values(end+1,1)=best; flags(end+1,1)=flag;
    iterations(end+1,1)=info.iterations; bestIndex=numel(values);
end
search=struct('initialTheta',starts,'theta',candidates,'negativeLogPosterior',values, ...
    'exitflag',flags,'iterations',iterations,'bestIndex',bestIndex, ...
    'globalOptimumGuaranteed',false);
[best,posterior]=marginal_objective(theta,data,prior,fs);
poles=(-2*pi*theta(1:2).*theta(3:4)+ ...
    1i*2*pi*theta(1:2).*sqrt(1-theta(3:4).^2)).';
residues=complex(zeros(2,2));
for j=1:2
    p=posterior(j); scale=data{j}.scale;
    residues(:,j)=(p.coefficientMean(1:2)+1i*p.coefficientMean(3:4))*scale;
    raw=2*real(exp(channels(j).relativeTime*poles.')*residues(:,j));
    centered=raw-mean(raw);
    channels(j).fittedRaw=raw;
    channels(j).fitted=centered;
    channels(j).offset=channels(j).removedMean-mean(raw);
    channels(j).fittedWithOffset=raw+channels(j).offset;
    channels(j).residual=channels(j).observed-centered;
    channels(j).fittedSpectrum=p.fittedSpectrum*scale;
    channels(j).relativeTimeError=norm(channels(j).residual)/norm(channels(j).observed);
    channels(j).relativeSpectrumError=norm(channels(j).spectrum-channels(j).fittedSpectrum) ...
        /max(norm(channels(j).spectrum),realmin);
    channels(j).noiseVarianceMean=p.b/(p.a-1)*scale^2;
    channels(j).conditionalResidueCovariance=p.coefficientCovariance*scale^2;
    channels(j).variancePosteriorShape=p.a;
    channels(j).variancePosteriorScale=p.b*scale^2;
end
shapes=residues.'; references=zeros(1,2);
for r=1:2
    [largest,references(r)]=max(abs(shapes(:,r)));
    if largest>0, shapes(:,r)=shapes(:,r)/shapes(references(r),r);
    else, shapes(:,r)=NaN; end
end
diagnostics=local_diagnostics(theta,lo,hi,objective,search);
diagnostics.modelOrder=2;
diagnostics.assumptions=['Two observable underdamped modes; force-free linear response; ', ...
    'constant modal parameters; independent white channel noise; common response type.'];
diagnostics.modelAdequacyVerified=false;
diagnostics.notes={ ...
    'Full positive spectrum is used. Extra modes, drift and colored noise can bias this two-mode fit.', ...
    'Conditional residue/noise moments hold poles fixed; local pole intervals are not global coverage guarantees.', ...
    'Common-response free decay supplies no absolute mass or stiffness scale.'};
if min((theta-lo)./span)<1e-6 || min((hi-theta)./span)<1e-6
    diagnostics.notes{end+1}='At least one estimate is close to its automatic prior support boundary.';
end
result=struct('frequencyHz',theta(1:2).','dampedFrequencyHz',abs(imag(poles))/(2*pi), ...
    'dampingRatio',theta(3:4).','dampingPercent',100*theta(3:4).', ...
    'poles',poles,'residues',residues,'modeShapes',shapes, ...
    'shapeReferenceChannel',references,'channels',channels,'prior',prior, ...
    'search',search,'diagnostics',diagnostics,'negativeLogPosterior',best, ...
    'samplingFrequency',fs,'startSampleTime',(first-1)/fs, ...
    'parameterOrder',{{'f1_Hz','f2_Hz','zeta1','zeta2'}}, ...
    'elapsedSeconds',toc(clockStart));
result.modalTable=table((1:2).',result.frequencyHz,result.dampedFrequencyHz, ...
    result.dampingRatio,result.dampingPercent, ...
    'VariableNames',{'Mode','NaturalFrequencyHz','DampedFrequencyHz','DampingRatio','DampingPercent'});
end

function index=grid_ceiling(value)
nearest=round(value);
if abs(value-nearest)<=16*eps(max(1,abs(value))), value=nearest; end
index=ceil(value);
end

function theta=decode(u,lo,hi)
u=max(-32,min(32,u));
theta=lo+(hi-lo)./(1+exp(-u));
[frequencies,order]=sort(theta(1:2));
theta=[frequencies,theta(2+order)];
end

function [value,out]=marginal_objective(theta,data,prior,fs)
value=1e100; out=struct([]);
if any(~isfinite(theta)) || any(theta<=prior.lower) || any(theta>=prior.upper) || theta(1)>=theta(2)
    return
end
lambda=-2*pi*theta(1:2).*theta(3:4)+1i*2*pi*theta(1:2).*sqrt(1-theta(3:4).^2);
value=0;
for j=1:2
    dc=data{j}; N=dc.N;
    plus=expm1(N*lambda/fs)./expm1(lambda/fs-1i*dc.omega)/sqrt(N);
    minus=expm1(N*conj(lambda)/fs)./expm1(conj(lambda)/fs-1i*dc.omega)/sqrt(N);
    P=[plus+minus,1i*(plus-minus)];
    D=sqrt(2)*[real(P);imag(P)];
    % Augmented QR avoids cancellation of Y'Y - Y'D*(D'D+ridge)^-1*D'Y.
    augmented=[D;eye(4)/prior.tau];
    [Q,R]=qr(augmented,0);
    w=R\(Q.'*[dc.Y;zeros(4,1)]);
    residual=dc.Y-D*w;
    b=prior.b0Normalized+.5*(sum(residual.^2)+sum(w.^2)/prior.tau^2);
    a=prior.a0+numel(dc.Y)/2;
    contribution=sum(log(abs(diag(R))))+4*log(prior.tau)+a*log(b);
    if ~isfinite(contribution), value=1e100; return; end
    value=value+contribution;
    if nargout>1
        invR=R\eye(4);
        out(j).coefficientMean=w;
        out(j).coefficientCovariance=(b/(a-1))*(invR*invR.');
        out(j).a=a; out(j).b=b;
        out(j).fittedSpectrum=P*w;
    end
end
end

function starts=automatic_starts(data,lo,hi,fs)
% Spectral peaks initialize only: all positive bins still enter inference.
candidateF=[]; scores=[]; dominant=zeros(1,2);
for j=1:2
    f=data{j}.frequencyHz; power=abs(data{j}.Z).^2;
    [peak,ix]=max(power); dominant(j)=f(ix);
    left=[-Inf;power(1:end-1)]; right=[power(2:end);-Inf];
    locations=find(power>=left & power>right);
    [~,order]=sort(power(locations),'descend'); locations=locations(order(1:min(4,end)));
    candidateF=[candidateF;f(locations)]; %#ok<AGROW>
    scores=[scores;power(locations)/max(peak,realmin)]; %#ok<AGROW>
end
[~,order]=sort(scores,'descend'); candidateF=candidateF(order);
resolution=fs/max([data{1}.N,data{2}.N]);
seeds=[];
for j=1:numel(candidateF)
    if isempty(seeds) || all(abs(seeds-candidateF(j))>.6*resolution)
        seeds(end+1)=candidateF(j); %#ok<AGROW>
    end
    if numel(seeds)==5, break; end
end
pairs=[];
for i=1:numel(seeds)
    for j=i+1:numel(seeds), pairs(end+1,:)=sort([seeds(i),seeds(j)]); end %#ok<AGROW>
end
% Seed unresolved close poles and several broader splits around each peak.
for center=unique(dominant)
    for gap=[.75*resolution,max(2*resolution,.03*center),.15*center]
        pairs(end+1,:)=[center-gap/2,center+gap/2]; %#ok<AGROW>
    end
end
pairs=[sort(dominant);pairs];
margin=max(resolution*.02,(hi(1)-lo(1))*1e-8);
pairs=max(lo(1)+margin,min(hi(1)-margin,pairs));
pairs=unique(pairs,'rows','stable');
pairs=pairs(diff(pairs,1,2)>margin,:);
if isempty(pairs), pairs=[fs/8,fs/4]; end
pairs=pairs(1:min(16,size(pairs,1)),:);
damping=[.01,.03;.05,.005;.003,.08;.04,.04];
starts=zeros(size(pairs,1)+8,4);
for j=1:size(pairs,1)
    starts(j,:)=[pairs(j,:),damping(mod(j-1,4)+1,:)];
end
for j=1:8
    starts(size(pairs,1)+j,:)=[pairs(mod(j-1,min(2,size(pairs,1)))+1,:),damping(mod(j-1,4)+1,:)];
end
starts=unique(starts,'rows','stable');
end

function diagnostics = local_diagnostics(theta,lo,hi,objective,search)
%LOCAL_DIAGNOSTICS Conservative numerical checks for local pole intervals.
% The density is defined in physical theta coordinates; scaling below only
% conditions derivatives. Passing these checks does not establish coverage.
theta=double(theta(:).');lo=double(lo(:).');hi=double(hi(:).');d=numel(theta);
diagnostics=struct('intervalsAccepted',false,'covariance',nan(d), ...
    'standardDeviation',nan(1,d),'ci95',nan(d,2),'status','not evaluated', ...
    'optimizerConverged',false,'scale',abs(theta),'scaledHessian',nan(d), ...
    'physicalHessian',nan(d),'scaledGradient',nan(d,1), ...
    'standardizedGradientNorm',NaN,'relativeHessianChanges',nan(1,2), ...
    'relativeSDChanges',nan(1,2),'scaledConditionNumber',NaN, ...
    'gapStandardDeviations',NaN,'separatedNearOptimalStarts',[], ...
    'finiteDifferenceSteps',nan(3,d),'candidateCI95',nan(d,2), ...
    'interpretation','Local Gaussian approximation conditional on the selected model, band and priors.');
issues={};
if d~=4 || any(~isfinite(theta)) || any(theta<=lo | theta>=hi) || theta(2)<=theta(1)
    diagnostics.status='Invalid or boundary pole estimate.';return
end
best=search.bestIndex;
diagnostics.optimizerConverged=isfinite(search.exitflag(best)) && search.exitflag(best)>0;
if ~diagnostics.optimizerConverged,issues{end+1}='Selected optimizer did not converge';end
scale=max(abs(theta),1e-8);diagnostics.scale=scale;
margin=min(theta-lo,hi-theta);
initialStep=min(2e-3*ones(1,d),margin./(4*scale));
initialStep(1:2)=min(initialStep(1:2),(theta(2)-theta(1))./(8*scale(1:2)));
steps=[initialStep;initialStep/2;initialStep/4];diagnostics.finiteDifferenceSteps=steps;
if any(steps(3,:)<1e-7)
    diagnostics.status=strjoin([issues,{'Insufficient room for reliable derivative stencils'}],'; ');return
end
H=nan(d,d,3);g=nan(d,3);C=nan(d,d,3);sd=nan(3,d);good=false(1,3);
try
    value=objective(theta);
    if ~isscalar(value) || ~isfinite(value),error('spectrum:Derivative','Nonfinite posterior value.');end
    for k=1:3
        [H(:,:,k),g(:,k)]=local_derivatives(theta,scale,steps(k,:),value,objective);
        eigval=eig(H(:,:,k));
        good(k)=all(isfinite(eigval)) && min(eigval)>1e-10*max(eigval) && min(eigval)>0;
        if good(k)
            R=chol(H(:,:,k));C(:,:,k)=R\(R'\eye(d));
            sd(k,:)=sqrt(diag(C(:,:,k))).'.*scale;
        end
    end
catch err
    diagnostics.status=strjoin([issues,{['Derivative evaluation failed: ',err.message]}],'; ');return
end
diagnostics.scaledHessian=H(:,:,3);
diagnostics.physicalHessian=H(:,:,3)./(scale.'*scale);
diagnostics.scaledGradient=g(:,3);
for k=1:2
    diagnostics.relativeHessianChanges(k)=norm(H(:,:,k+1)-H(:,:,k),'fro')/max(norm(H(:,:,k+1),'fro'),realmin);
    if good(k) && good(k+1)
        diagnostics.relativeSDChanges(k)=max(abs(sd(k+1,:)-sd(k,:))./max(sd(k+1,:),realmin));
    end
end
if ~all(good)
    issues{end+1}='Local Hessian is not positive definite and well conditioned at every step';
else
    diagnostics.scaledConditionNumber=cond(H(:,:,3));
    diagnostics.standardizedGradientNorm=sqrt(max(0,g(:,3)'*C(:,:,3)*g(:,3)));
    if diagnostics.standardizedGradientNorm>0.05
        issues{end+1}='Posterior gradient is too large relative to local uncertainty';
    end
    if any(diagnostics.relativeHessianChanges>0.05) || any(diagnostics.relativeSDChanges>0.10)
        issues{end+1}='Local curvature or uncertainty is unstable across derivative steps';
    end
    covariance=C(:,:,3).*(scale.'*scale);covariance=(covariance+covariance')/2;
    standardDeviation=sqrt(diag(covariance)).';
    ci=[theta.'-1.95996398454005*standardDeviation.',theta.'+1.95996398454005*standardDeviation.'];
    diagnostics.candidateCI95=ci;
    if any(ci(:,1)<=lo.' | ci(:,2)>=hi.')
        issues{end+1}='Local intervals cross the pole-prior support';
    end
    gapVariance=covariance(1,1)+covariance(2,2)-2*covariance(1,2);
    diagnostics.gapStandardDeviations=(theta(2)-theta(1))/sqrt(max(gapVariance,realmin));
    if diagnostics.gapStandardDeviations<=1.95996398454005
        issues{end+1}='Local uncertainty crosses the ordered-frequency boundary';
    end
    values=search.negativeLogPosterior(:);
    candidates=find(isfinite(values) & values<=value+2);
    competing=[];
    for j=candidates.'
        alternative=search.theta(j,:);
        if all(isfinite(alternative))
            delta=(alternative-theta)./scale;
            distanceSquared=delta*H(:,:,3)*delta';
            if distanceSquared>16,competing(end+1)=j;end %#ok<AGROW>
        end
    end
    diagnostics.separatedNearOptimalStarts=competing;
    if ~isempty(competing),issues{end+1}='Separated near-optimal multistart candidates remain';end
    if isempty(issues)
        diagnostics.intervalsAccepted=true;diagnostics.covariance=covariance;
        diagnostics.standardDeviation=standardDeviation;diagnostics.ci95=ci;
    end
end
if isempty(issues)
    diagnostics.status='Local numerical checks passed; Gaussian intervals are approximate.';
else
    diagnostics.status=strjoin(issues,'; ');
end
end

function [H,g]=local_derivatives(theta,scale,h,value,objective)
d=numel(theta);H=zeros(d);g=zeros(d,1);
for i=1:d
    ei=zeros(1,d);ei(i)=h(i)*scale(i);
    plus=objective(theta+ei);minus=objective(theta-ei);
    g(i)=(plus-minus)/(2*h(i));H(i,i)=(plus-2*value+minus)/h(i)^2;
    for j=1:i-1
        ej=zeros(1,d);ej(j)=h(j)*scale(j);
        H(i,j)=(objective(theta+ei+ej)-objective(theta+ei-ej) ...
            -objective(theta-ei+ej)+objective(theta-ei-ej))/(4*h(i)*h(j));
        H(j,i)=H(i,j);
    end
end
if any(~isfinite(H(:))) || any(~isfinite(g))
    error('spectrum:Derivative','Nonfinite finite-difference derivative.');
end
end

