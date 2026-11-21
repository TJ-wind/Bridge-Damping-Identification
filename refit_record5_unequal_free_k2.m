function physicalFit = refit_record5_unequal_free_k2(bridgeResponse,exciterResponse,fs, ...
    startTimes,m2,hasInerter,responseOrders,inputScaleToSI,frequencyBandHz, ...
    priorLower,priorUpper,tau,a0,b0,statePriorScale,initialParameters,initialState, ...
    searchSeeds,populationSize,maxIterations,maxStallIterations,localStarts,localMaxIterations)
%REFIT_RECORD5_UNEQUAL_FREE_K2 两路已选段响应的独立物理参数拟合。
% 保留旧函数名称；已去除 test5、旧结果 MAT 和固定窗口依赖。
% 本函数不调用或修改 main_function，不修改其 result，不自动保存文件。
%
% 输入均为明确变量；可在 main_physical_refit.m 的输入区逐项设置。
% bridgeResponse,exciterResponse: 已选段、尚未去均值的两路响应向量。
% fs: 共同采样频率 Hz；startTimes: 两路首个实际采样点的绝对时间 s。
% m2: 已知振子质量 kg；hasInerter=false 时 beta 严格固定为 0。
% responseOrders: [桥梁阶次,振子阶次]，0=位移、1=速度、2=加速度。
% inputScaleToSI: 两路乘数，分别转换到 m、m/s 或 m/s^2，可含已知极性。
% frequencyBandHz: [下限,上限]，每路保留自身 N 对应的正频率 DFT 频点。
% priorLower,priorUpper,initialParameters 的列顺序均为：
%   [rho, fE, zetaB, zetaE, beta, rN, fB]
%   E=m2/(1-beta),rho=mB/E,beta=b/E,rN=sigmaD/sigmaB;
%   fB=sqrt(kB/mB)/(2*pi),fE=sqrt(kE/E)/(2*pi).
%   fB,fE,zetaB,zetaE 是元件参数，不能直接称为耦合模态频率/阻尼。
% 已知调谐参数：令第2列 fE 和第5列 beta 的上下界相等即可固定。
% 固定 beta 给定 b；同时固定 fE 给定 kE=(m2+b)*(2*pi*fE)^2。
% initialParameters 可有多行；initialState=[] 或 [xB0,xD0,vB0,vD0]。
% 状态共同参考时刻为 min(startTimes)，不同起点精确传播至各路首点。
% statePriorScale: 4个正缩放系数，u=diag(scale)*[xB0;xD0;vB0;vD0]。
% 原速度先验为 [2*pi*fRef,2*pi*fRef,1,1]；加速度为 [wRef^2,wRef^2,wRef,wRef]。
% 混合响应时，缩放系数应将 u 转为桥梁响应的 SI 单位；
% rN 此时是两路 SI 噪声标准差之比，单位随两路响应类型确定。
% u|sigmaB^2~N(0,tau^2*sigmaB^2*I), sigmaB^2~IG(a0,b0).
% 机械参数在上述物理坐标中使用有界均匀先验；优化变换不加入 Jacobian。
% 四个初始状态保持自由，经增广 QR 剖面求最优；仅噪声方差被积分。
% 目标与旧版相同：L2*log(rN)+(a0+(L1+L2+4)/2)*log(b0+Qmin/2)。
% 搜索：联合状态局部、多随机种子 PSO/DE、多候选 SQP、最终联合阶段。
% 输出 physicalFit 与原谱识别 result 分开；正参数约束不等于反演已验证。

if nargin~=23
    error('physicalRefit:ExplicitInputs', ...
        '需要23个明确输入参数。请运行 main_physical_refit.m，或查看 help refit_record5_unequal_free_k2。');
end
validateattributes(bridgeResponse,{'numeric'},{'real','vector','finite','nonempty'});
validateattributes(exciterResponse,{'numeric'},{'real','vector','finite','nonempty'});
validateattributes(fs,{'numeric'},{'real','scalar','positive','finite'});
validateattributes(startTimes,{'numeric'},{'real','vector','numel',2,'finite'});
validateattributes(m2,{'numeric'},{'real','scalar','positive','finite'});
assert(isscalar(hasInerter) && (islogical(hasInerter) || ismember(hasInerter,[0,1])), ...
    'physicalRefit:Topology','hasInerter must be true or false.');
validateattributes(responseOrders,{'numeric'},{'real','vector','numel',2,'integer','>=',0,'<=',2});
validateattributes(inputScaleToSI,{'numeric'},{'real','vector','numel',2,'finite','nonzero'});
validateattributes(frequencyBandHz,{'numeric'},{'real','vector','numel',2,'finite','nonnegative'});
assert(frequencyBandHz(2)>frequencyBandHz(1) && frequencyBandHz(2)<=fs/2, ...
    'physicalRefit:FrequencyBand','Require 0<=lower<upper<=fs/2.');
validateattributes(priorLower,{'numeric'},{'real','vector','numel',7,'finite'});
validateattributes(priorUpper,{'numeric'},{'real','vector','numel',7,'finite'});
priorLower=double(priorLower(:).');priorUpper=double(priorUpper(:).');
positive=[1:4,6,7];
freePositive=[1,3,4,6,7];
assert(all(priorLower(positive)>0) && all(priorUpper(freePositive)>priorLower(freePositive)) ...
    && priorUpper(2)>=priorLower(2), ...
    'physicalRefit:PriorBounds','Positive coordinates require 0<lower<upper; known fE may use equal positive bounds.');
if hasInerter
    assert(priorLower(5)>=0 && priorUpper(5)<1 && priorUpper(5)>=priorLower(5), ...
        'physicalRefit:BetaBounds','With an inerter require 0<=betaLower<=betaUpper<1; equal bounds fix known inertance.');
else
    assert(priorLower(5)==0 && priorUpper(5)==0, ...
        'physicalRefit:NoInerterBounds','Without an inerter set priorLower(5)=priorUpper(5)=0.');
end
for v={tau,a0,b0}
    validateattributes(v{1},{'numeric'},{'real','scalar','positive','finite'});
end
validateattributes(statePriorScale,{'numeric'},{'real','vector','numel',4,'positive','finite'});
validateattributes(initialParameters,{'numeric'},{'real','2d','ncols',7,'finite','nonempty'});
initialParameters=double(initialParameters);
assert(all(initialParameters>=priorLower & initialParameters<=priorUpper,'all'), ...
    'physicalRefit:InitialParameters','Every initial parameter row must lie in the explicit prior bounds.');
if ~isempty(initialState)
    validateattributes(initialState,{'numeric'},{'real','vector','numel',4,'finite'});
end
validateattributes(searchSeeds,{'numeric'},{'real','vector','integer','nonnegative','<=',2^32-1,'finite','nonempty'});
for v={populationSize,maxIterations,maxStallIterations,localStarts,localMaxIterations}
    validateattributes(v{1},{'numeric'},{'real','scalar','integer','positive','finite'});
end
assert(populationSize>=max(4,size(initialParameters,1)+1), ...
    'physicalRefit:Population','Population must include all initial rows, the joint candidate, and at least four members.');
assert(license('test','Optimization_Toolbox') && exist('fmincon','file')==2, ...
    'physicalRefit:OptimizationToolbox','The original SQP/joint stages require Optimization Toolbox (fmincon).');
usePSO=license('test','GADS_Toolbox') && exist('particleswarm','file')==2;
savedRng=rng;restoreRng=onCleanup(@()rng(savedRng));
startTimes=double(startTimes(:).');responseOrders=double(responseOrders(:).');
inputScaleToSI=double(inputScaleToSI(:).');statePriorScale=double(statePriorScale(:));
referenceStartTime=min(startTimes);
signals={double(bridgeResponse(:))*inputScaleToSI(1),double(exciterResponse(:))*inputScaleToSI(2)};
data=struct('fs',fs,'referenceStartTime',referenceStartTime, ...
    'startOffsets',startTimes-referenceStartTime,'responseOrders',responseOrders,'channels',{{}});
timeChannels=cell(1,2);timeAbsoluteChannels=cell(1,2);yChannels=cell(1,2);
removedMeans=zeros(1,2);N=zeros(1,2);
for c=1:2
    N(c)=numel(signals{c});
    assert(N(c)>=8,'physicalRefit:SampleCount','Each selected channel needs at least eight samples.');
    removedMeans(c)=mean(signals{c});yChannels{c}=signals{c}-removedMeans(c);
    assert(norm(yChannels{c})>0,'physicalRefit:ConstantSignal','Selected channels must contain a varying response.');
    timeChannels{c}=(0:N(c)-1).'/fs;
    timeAbsoluteChannels{c}=startTimes(c)+timeChannels{c};
    k=(2:floor((N(c)-1)/2)+1).';frequencies=(k-1)*fs/N(c);
    retain=frequencies>=frequencyBandHz(1) & frequencies<=frequencyBandHz(2);
    k=k(retain);frequencies=frequencies(retain);
    assert(~isempty(k),'physicalRefit:EmptyBand','The selected band contains no positive non-Nyquist bins for channel %d.',c);
    F=fft(yChannels{c})/sqrt(N(c));
    data.channels{c}=struct('N',N(c),'dt',1/fs,'k',k,'frequencies',frequencies, ...
        'Y',sqrt(2)*[real(F(k));imag(F(k))]);
end
active=find(priorUpper>priorLower);d=numel(active);
% Original transforms: rho,fE,zetaB,zetaE,rN logarithmic; beta,fB linear.
logMask=ismember(active,[1:4,6]);
lowerQ=coordinate(priorLower(active),logMask);upperQ=coordinate(priorUpper(active),logMask);
encode=@(p)(coordinate(p(active),logMask)-lowerQ)./(upperQ-lowerQ);
decode=@(q)expand(uncoordinate(lowerQ+min(max(q(:).',0),1).*(upperQ-lowerQ),logMask), ...
    active,priorLower,priorUpper);
evaluationCount=0;timer=tic;profile=@(q)evaluate(decode(q));
localOptions=optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxIterations',localMaxIterations,'MaxFunctionEvaluations',max(14000,20*localMaxIterations), ...
    'StepTolerance',1e-10,'OptimalityTolerance',1e-7,'FunctionTolerance',1e-10, ...
    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-6);
jointOptions=optimoptions(localOptions,'SpecifyObjectiveGradient',true);
candidateQ=zeros(size(initialParameters,1),d);candidateValue=zeros(size(initialParameters,1),1);
for j=1:size(initialParameters,1)
    candidateQ(j,:)=encode(initialParameters(j,:));candidateValue(j)=profile(candidateQ(j,:));
end
[baselineValue,first]=min(candidateValue);
assert(baselineValue<1e90,'physicalRefit:InitialModel', ...
    'None of the initial parameter rows gives two stable, non-aliased underdamped modes.');
baselineQ=candidateQ(first,:);[~,baselineOut]=profile(baselineQ);
if isempty(initialState),initialU=baselineOut.u0;else,initialU=initialState(:).*statePriorScale;end
% Numerical optimizer scaling only; the supplied state prior is unchanged.
jointScale=max(abs(initialU),max(sqrt(mean(yChannels{1}.^2)),realmin)*ones(4,1));
jointStages=cell(1,2);jointStages{1}=run_joint(baselineQ,initialU,'initial_joint_stage');
candidateQ(end+1,:)=jointStages{1}.profileQ;candidateValue(end+1,1)=jointStages{1}.profileValue;
globalRuns=cell(1,numel(searchSeeds));
for run=1:numel(searchSeeds)
    rng(searchSeeds(run),'twister');
    fixedRows=[candidateQ(1:size(initialParameters,1),:);jointStages{1}.profileQ];
    population=[fixedRows;rand(populationSize-size(fixedRows,1),d)];
    fprintf('Physical refit: global start %d/%d, seed %d.\n',run,numel(searchSeeds),searchSeeds(run));
    if usePSO
        swarmOptions=optimoptions('particleswarm','Display','off','SwarmSize',populationSize, ...
            'MaxIterations',maxIterations,'MaxStallIterations',maxStallIterations, ...
            'FunctionTolerance',1e-8,'InitialSwarmMatrix',population, ...
            'UseParallel',false,'UseVectorized',false);
        [q,value,flag,output,points]=particleswarm(profile,d,zeros(1,d),ones(1,d),swarmOptions);
        algorithm='particleswarm';
    else
        [q,value,flag,output,points]=physical_refit_support.differential_evolution( ...
            profile,population,maxIterations,maxStallIterations);
        algorithm='differential_evolution_rand1bin';
    end
    globalRuns{run}=struct('algorithm',algorithm,'seed',searchSeeds(run), ...
        'parameters',decode(q),'nll',value,'exitflag',flag,'output',output);
    candidateQ=[candidateQ;q;points.X];candidateValue=[candidateValue;value;points.Fval(:)]; %#ok<AGROW>
end
finite=isfinite(candidateValue) & candidateValue<1e90;
candidateQ=candidateQ(finite,:);candidateValue=candidateValue(finite);
[~,order]=sort(candidateValue);selected=[];
for idx=order(:).'
    if isempty(selected) || all(vecnorm(candidateQ(selected,:)-candidateQ(idx,:),2,2)>1e-3)
        selected(end+1)=idx; %#ok<AGROW>
    end
    if numel(selected)==localStarts,break;end
end
if numel(selected)<localStarts
    for idx=order(:).'
        if ~ismember(idx,selected),selected(end+1)=idx;end %#ok<AGROW>
        if numel(selected)==localStarts,break;end
    end
end
localRuns=cell(1,numel(selected));
for j=1:numel(selected)
    q0=candidateQ(selected(j),:);
    record=struct('initialValue',candidateValue(selected(j)),'parameters',decode(q0), ...
        'nll',candidateValue(selected(j)),'exitflag',-99,'output',struct(),'error','');
    try
        [q,value,flag,output]=fmincon(profile,q0,[],[],[],[],zeros(1,d),ones(1,d),[],localOptions);
        record.parameters=decode(q);record.nll=value;record.exitflag=flag;record.output=output;
        candidateQ(end+1,:)=q;candidateValue(end+1,1)=value; %#ok<AGROW>
    catch err
        record.error=sprintf('%s: %s',err.identifier,err.message);
    end
    localRuns{j}=record;
end
[~,idx]=min(candidateValue);bestQ=candidateQ(idx,:);[~,bestOut]=profile(bestQ);
finalU=bestOut.u0;if ~isempty(initialState),finalU=initialState(:).*statePriorScale;end
jointStages{2}=run_joint(bestQ,finalU,'final_joint_stage');
candidateQ(end+1,:)=jointStages{2}.profileQ;candidateValue(end+1,1)=jointStages{2}.profileValue;
[~,idx]=min(candidateValue);bestQ=candidateQ(idx,:);p=decode(bestQ);
[bestValue,out]=evaluate(p);
assert(bestValue<=baselineValue+1e-8,'physicalRefit:Regression','The retained objective regressed from its initial value.');
m=out.model;w=out.stateMap\out.u0;cx=w(1:2)+1i*w(3:4);
yhatChannels=cell(1,2);centeredYhatChannels=cell(1,2);fittedWithOffsetsChannels=cell(1,2);
residualChannels=cell(1,2);dataFourier=cell(1,2);modelFourier=cell(1,2);
complexResidual=cell(1,2);responseResidues=complex(zeros(2,2));
timeNRMSE=zeros(1,2);directTimeNRMSE=zeros(1,2);selectedNRMSE=zeros(1,2);
fourierVerificationError=zeros(1,2);
for c=1:2
    ch=data.channels{c};
    responseResidues(:,c)=cx.*m.displacementShapes(c,:).'.*(m.lambda.^responseOrders(c)).* ...
        exp(m.lambda*data.startOffsets(c));
    yhatChannels{c}=2*real(exp(timeChannels{c}*m.lambda.')*responseResidues(:,c));
    centeredYhatChannels{c}=yhatChannels{c}-mean(yhatChannels{c});
    fittedWithOffsetsChannels{c}=centeredYhatChannels{c}+removedMeans(c);
    residualChannels{c}=signals{c}-fittedWithOffsetsChannels{c};
    F=fft(yChannels{c})/sqrt(ch.N);G=fft(yhatChannels{c})/sqrt(ch.N);
    dataFourier{c}=F(ch.k);modelFourier{c}=G(ch.k);complexResidual{c}=F(ch.k)-G(ch.k);
    timeNRMSE(c)=norm(residualChannels{c})/norm(yChannels{c});
    directTimeNRMSE(c)=norm(yChannels{c}-yhatChannels{c})/norm(yChannels{c});
    selectedNRMSE(c)=norm(complexResidual{c})/max(norm(F(ch.k)),realmin);
    exactStack=sqrt(2)*[real(G(ch.k));imag(G(ch.k))];
    designed=out.Xraw(out.rowIndices{c},:)*out.u0;
    fourierVerificationError(c)=norm(exactStack-designed)/max(norm(exactStack),realmin);
end
assert(all(fourierVerificationError<1e-8),'physicalRefit:FourierVerification','Time response and finite DFT design disagree.');
lowerHit=false(1,7);upperHit=false(1,7);
lowerHit(active)=bestQ<=1e-6;upperHit(active)=bestQ>=1-1e-6;
fixedCoordinates=priorLower==priorUpper;
parameterNames={'rho','fE_Hz','zetaB','zetaE','beta','rN','fB_Hz'};
parameterTable=table(string(parameterNames(:)),p(:),priorLower(:),priorUpper(:), ...
    fixedCoordinates(:),lowerHit(:),upperHit(:), ...
    'VariableNames',{'Parameter','Estimate','Lower','Upper','Fixed','AtLowerBound','AtUpperBound'});
mechanicalNames={'mB_kg','m2_kg','b_kg','kB_N_per_m','kE_N_per_m','cB_Ns_per_m','cE_Ns_per_m'};
mechanicalValues=[m.m1,m.m2,m.b,m.k1,m.k2,m.c1,m.c2];
physicalParameterTable=table(string(mechanicalNames(:)),mechanicalValues(:), ...
    [false;true;fixedCoordinates(5);false;fixedCoordinates(2)&&fixedCoordinates(5);false;false], ...
    'VariableNames',{'Parameter','Estimate','Fixed'});
[~,massFlag]=chol(m.M);[~,stiffnessFlag]=chol(m.K);[~,dampingFlag]=chol(m.C);
positiveMechanical=all(mechanicalValues([1,2,4:7])>0) && m.b>=0 && ...
    massFlag==0 && stiffnessFlag==0 && dampingFlag==0;
anyLocalConverged=any(cellfun(@(r)r.exitflag>0,localRuns)) || ...
    any(cellfun(@(r)r.exitflag>0,jointStages));
convergedAtBest=any(cellfun(@(r)r.exitflag>0 && abs(r.nll-bestValue)<=1e-7*max(1,abs(bestValue)),localRuns)) || ...
    any(cellfun(@(r)r.exitflag>0 && abs(r.profileValue-bestValue)<=1e-7*max(1,abs(bestValue)),jointStages));
if any(lowerHit|upperHit)
    status='bounded_candidate';
elseif convergedAtBest
    status='local_stationary_candidate';
else
    status='search_limit_candidate';
end
physicalFit=struct('status',status,'fullParameters',p,'parameters',p(1:6), ...
    'fullParameterNames',{parameterNames},'model',m,'lambda',m.lambda,'coupled_theta',m.theta, ...
    'frequencyHz',m.theta(1:2),'dampingRatio',m.theta(3:4),'dampingPercent',100*m.theta(3:4), ...
    'modeShapes',m.displacementShapes,'u0',out.u0,'state0',out.state0, ...
    'stateOrder',{{'xB0','xD0','vB0','vD0'}},'statePriorScale',statePriorScale, ...
    'referenceStartTime',referenceStartTime,'startTimes',startTimes,'responseOrders',responseOrders, ...
    'inputScaleToSI',inputScaleToSI,'fs',fs,'hasInerter',logical(hasInerter),'fixedExciterMass',m2, ...
    'frequencyBandHz',frequencyBandHz,'knownInputs',struct('m2',m2,'hasInerter',logical(hasInerter)), ...
    'prior',struct('lower',priorLower,'upper',priorUpper,'tau',tau,'a0',a0,'b0',b0), ...
    'data',data,'observedSIChannels',{signals},'yChannels',{yChannels},'removedMeans',removedMeans, ...
    'timeChannels',{timeChannels},'timeAbsoluteChannels',{timeAbsoluteChannels}, ...
    'yhatChannels',{yhatChannels},'centeredYhatChannels',{centeredYhatChannels}, ...
    'fittedWithOffsetsChannels',{fittedWithOffsetsChannels},'residualChannels',{residualChannels}, ...
    'dataFourier',{dataFourier},'modelFourier',{modelFourier},'complexResidual',{complexResidual}, ...
    'displacementModalCoefficients',cx,'responseResidues',responseResidues, ...
    'centeredTimeNRMSE',timeNRMSE,'timeNRMSE',timeNRMSE,'directTimeNRMSE',directTimeNRMSE, ...
    'selectedNRMSE',selectedNRMSE,'nll',bestValue,'baselineNll',baselineValue, ...
    'noiseRatio',p(6),'noiseVariancePosteriorShape',out.a,'noiseVariancePosteriorScale',out.b, ...
    'noiseSDPosteriorMeanVarianceSqrt',sqrt(out.b/(out.a-1))*[1,p(6)], ...
    'parameterTable',parameterTable,'physicalParameterTable',physicalParameterTable, ...
    'boundHits',struct('lower',lowerHit,'upper',upperHit,'fixed',fixedCoordinates), ...
    'mechanicalPositive',positiveMechanical,'physicalParametersValidated',false, ...
    'credibleIntervalsAccepted',false,'globalOptimalityCertified',false, ...
    'convergence',struct('anyLocalConverged',anyLocalConverged,'convergedAtBest',convergedAtBest, ...
        'globalRuns',{globalRuns},'localRuns',{localRuns},'jointStages',{jointStages}), ...
    'searchSettings',struct('seeds',searchSeeds,'populationSize',populationSize, ...
        'maxIterations',maxIterations,'maxStallIterations',maxStallIterations, ...
        'localStarts',localStarts,'localMaxIterations',localMaxIterations, ...
        'initialParameters',initialParameters,'initialState',initialState,'usedPSO',usePSO), ...
    'diagnostics',struct('fourierVerificationError',fourierVerificationError, ...
        'modeEquationRelativeError',m.modeEquationRelativeError,'evaluationCount',evaluationCount, ...
        'initialStatePriorPenalty',out.priorPenalty,'weightedResidualRSS',out.residualRSS), ...
    'elapsedSeconds',toc(timer), ...
    'objectiveDefinition','L2*log(rN)+(a0+(L1+L2+4)/2)*log(b0+Qmin/2); four shared states profiled, noise variance integrated; no design determinant or coordinate Jacobian.', ...
    'priorDefinition','Uniform bounded density in rho,fE,zetaB,zetaE,beta,rN,fB; u|sigmaB^2 is isotropic Gaussian with the explicitly supplied state scale.', ...
    'statePriorTimeMeaning','The state prior is imposed at min(startTimes); changing this reference changes the stated prior.', ...
    'primaryTimeComparison','Measured and fitted finite-record means are projected separately for each selected channel.', ...
    'validationMeaning','Positive mechanics satisfy the chosen passive topology. Search convergence and parameter positivity do not establish uniqueness or calibrated physical uncertainty.');
physicalFit.diagnostics.status=status;
if fixedCoordinates(5)
    physicalFit.knownInputs.inertanceKg=m.b;
end
if fixedCoordinates(2)
    physicalFit.knownInputs.fixedBaseExciterFrequencyHz=p(2);
end
if fixedCoordinates(2) && fixedCoordinates(5)
    physicalFit.knownInputs.attachmentStiffnessNPerM=m.k2;
end
physicalFit.diagnostics.mechanicalPositive=positiveMechanical;
physicalFit.diagnostics.boundHits=physicalFit.boundHits;
physicalFit.diagnostics.anyLocalConverged=anyLocalConverged;
physicalFit.diagnostics.convergedAtBest=convergedAtBest;
physicalFit.diagnostics.physicalParametersValidated=false;
physicalFit.diagnostics.globalOptimalityCertified=false;
fprintf('Physical refit: %s; time NRMSE %.4f%% / %.4f%%; selected-DFT NRMSE %.4f%% / %.4f%%.\n', ...
    status,100*timeNRMSE,100*selectedNRMSE);

    function [value,details]=evaluate(theta,u)
        evaluationCount=evaluationCount+1;
        if nargin<2
            [value,details]=physical_refit_support.objective(theta,data,m2,tau,a0,b0,statePriorScale);
        else
            [value,details]=physical_refit_support.objective(theta,data,m2,tau,a0,b0,statePriorScale,u);
        end
    end
    function stage=run_joint(q0,u0,label)
        z0=[q0,(u0./jointScale).'];initialValue=joint_objective(z0);
        stage=struct('label',label,'initialValue',initialValue,'nll',initialValue, ...
            'exitflag',-99,'output',struct(),'error','', ...
            'profileQ',q0,'profileValue',profile(q0),'u0',u0);
        try
            [z,value,flag,output]=fmincon(@joint_objective,z0,[],[],[],[], ...
                [zeros(1,d),-inf(1,4)],[ones(1,d),inf(1,4)],[],jointOptions);
            stage.nll=value;stage.exitflag=flag;stage.output=output;
            stage.profileQ=z(1:d);stage.profileValue=profile(z(1:d));
            stage.u0=z(d+1:end).'.*jointScale;
        catch err
            stage.error=sprintf('%s: %s',err.identifier,err.message);
        end
    end
    function [value,gradient]=joint_objective(z)
        z=z(:).';q=z(1:d);u=z(d+1:end).'.*jointScale;
        [value,details]=evaluate(decode(q),u);
        if nargout>1
            gradient=zeros(d+4,1);
            for jj=1:d
                qp=q;qm=q;qp(jj)=min(1,q(jj)+1e-6);qm(jj)=max(0,q(jj)-1e-6);
                gradient(jj)=(evaluate(decode(qp),u)-evaluate(decode(qm),u))/(qp(jj)-qm(jj));
            end
            if isfield(details,'gradientU'),gradient(d+1:end)=details.gradientU(:).*jointScale;end
        end
    end
end

function p=expand(values,active,lower,upper)
p=lower;p(active)=min(max(values,lower(active)),upper(active));
end
function x=coordinate(x,mask)
x(mask)=log(x(mask));
end
function x=uncoordinate(x,mask)
x(mask)=exp(x(mask));
end
