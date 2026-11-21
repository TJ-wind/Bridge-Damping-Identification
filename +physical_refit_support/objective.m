function [value,out] = objective(p,data,m2,tau,a0,b0,statePriorScale,uOptional)
%OBJECTIVE Original integrated-noise / profiled shared-state objective.
% u=diag(statePriorScale)*[xB0;xD0;vB0;vD0] at data.referenceStartTime.
% u|sigmaB^2 ~ N(0,tau^2*sigmaB^2*I4), sigmaB^2 ~ IG(a0,b0).
% L2*log(rN)+(a0+(L1+L2+4)/2)*log(b0+Q/2), no design determinant.
% Distinct window starts enter exp(lambda*offset); each response derivative
% enters lambda^responseOrder. Same-start velocity is the original target.
value=1e100; out=struct();
try
    m=physical_refit_support.model(p,m2);
    if ~m.twoUnderdampedPairs,return;end
    lambda=m.lambda; Phi=m.displacementShapes;
    if any(abs(imag(lambda))/(2*pi)>=data.fs/2),return;end
    B=2*[real(Phi),-imag(Phi); ...
        real(Phi.*lambda.'),-imag(Phi.*lambda.')];
    scaling=statePriorScale(:);
    T=diag(scaling)*B;
    if rcond(T)<1e-14,return;end
    blocks=cell(1,2); observations=cell(1,2); L=zeros(1,2);
    for c=1:2
        ch=data.channels{c};
        D=finite_dft_design(lambda,ch.k,ch.N,1/data.fs);
        shape=Phi(c,:).*(lambda.'.^data.responseOrders(c)).* ...
            exp(lambda.'*data.startOffsets(c));
        map=[diag(real(shape)),-diag(imag(shape)); ...
            diag(imag(shape)),diag(real(shape))];
        blocks{c}=D*map/T;
        observations{c}=ch.Y;
        L(c)=numel(ch.Y);
    end
    Xraw=vertcat(blocks{:}); Yraw=vertcat(observations{:});
    rows={1:L(1),L(1)+(1:L(2))};
    X=Xraw; Y=Yraw; rN=p(6);
    X(rows{2},:)=X(rows{2},:)/rN;
    Y(rows{2})=Y(rows{2})/rN;
    if nargin<8 || isempty(uOptional)
        [Qmat,R]=qr([X;eye(4)/tau],0);
        u0=R\(Qmat.'*[Y;zeros(4,1)]);
        profiled=true;
    else
        u0=uOptional(:);
        assert(numel(u0)==4 && all(isfinite(u0)));
        profiled=false;
    end
    residual=Y-X*u0;
    penalty=sum(u0.^2)/tau^2;
    Q=sum(residual.^2)+penalty;
    a=a0+(numel(Y)+4)/2; b=b0+Q/2;
    value=L(2)*log(rN)+a*log(b);
    if ~isfinite(value),value=1e100;return;end
    if nargout>1
        out=struct('X',X,'Y',Y,'Xraw',Xraw,'Yraw',Yraw,'u0',u0, ...
            'state0',u0./scaling,'stateScaling',scaling, ...
            'stateMap',T,'physicalStateMap',B,'model',m, ...
            'residualRSS',sum(residual.^2),'priorPenalty',penalty, ...
            'Q',Q,'a',a,'b',b,'L',L,'rowIndices',{rows}, ...
            'noiseRatio',rN,'noiseLogDetCorrection',L(2)*log(rN), ...
            'gradientU',(a/b)*(-X.'*residual+u0/tau^2), ...
            'isProfiled',profiled);
    end
catch err
    value=1e100;
    if nargout>1
        out=struct('failureIdentifier',err.identifier,'failureMessage',err.message);
    end
end
end

function D=finite_dft_design(lambda,k,N,dt)
% fft(y)/sqrt(N), stacked sqrt(2)*[real;imag], positive non-Nyquist bins.
k=k(:); P=complex(zeros(numel(k),4));
for j=1:2
    gp=expm1(N*lambda(j)*dt)./expm1(lambda(j)*dt-1i*2*pi*(k-1)/N)/sqrt(N);
    gm=expm1(N*conj(lambda(j))*dt)./expm1(conj(lambda(j))*dt-1i*2*pi*(k-1)/N)/sqrt(N);
    P(:,j)=gp+gm;
    P(:,j+2)=1i*(gp-gm);
end
D=sqrt(2)*[real(P);imag(P)];
end
