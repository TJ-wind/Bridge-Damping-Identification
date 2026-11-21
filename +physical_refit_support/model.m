function model = model(p,m2)
%MODEL Original passive bridge--inerter--exciter mechanical parameterization.
% p = [rho,fE,zetaB,zetaE,beta,rN,fB], E=m2+b=m2/(1-beta).
% rho=mB/E; fB=sqrt(kB/mB)/(2*pi); fE=sqrt(kE/E)/(2*pi).
% These two element frequencies/damping ratios are not coupled modal values.
validateattributes(p,{'numeric'},{'real','vector','numel',7,'finite'});
validateattributes(m2,{'numeric'},{'real','scalar','positive','finite'});
p=double(p(:).');
assert(all(p([1:4,6,7])>0) && p(5)>=0 && p(5)<1, ...
    'physicalRefit:Support','Require positive rho,fE,zetaB,zetaE,rN,fB and 0<=beta<1.');
r=p(1); fE=p(2); zB=p(3); zE=p(4); beta=p(5);
E=m2/(1-beta); wB=2*pi*p(7); wE=2*pi*fE;
Mn=[r+beta,-beta;-beta,1];
Cn=[2*zB*r*wB+2*zE*wE,-2*zE*wE;-2*zE*wE,2*zE*wE];
Kn=[r*wB^2+wE^2,-wE^2;-wE^2,wE^2];
M=E*Mn; C=E*Cn; K=E*Kn;
assert(all(isfinite([M(:);C(:);K(:)])),'physicalRefit:Overflow','Physical matrices overflowed.');
A=[zeros(2),eye(2);-Mn\Kn,-Mn\Cn];
[V,L]=eig(A); allPoles=diag(L);
selected=find(imag(allPoles)<-1e-10*max(1,norm(A,2)));
valid=numel(selected)==2 && all(real(allPoles)<0);
model=struct('parameters',p,'known',struct('f1',p(7),'E',E), ...
    'M',M,'C',C,'K',K,'A',A,'E',E, ...
    'm1',r*E,'m2',(1-beta)*E,'b',beta*E, ...
    'k1',r*E*wB^2,'k2',E*wE^2,'c1',2*zB*r*E*wB,'c2',2*zE*E*wE, ...
    'uncoupledBridgeFrequencyHz',p(7),'fixedBaseExciterFrequencyHz',fE, ...
    'uncoupledBridgeDamping',zB,'fixedBaseExciterDamping',zE, ...
    'allPoles',allPoles,'twoUnderdampedPairs',valid,'lambda',[],'theta',[], ...
    'displacementShapes',[],'shapeReference',[],'shapeRatios',[], ...
    'massRatioMeaning','rho=mB/(m2+b); not a prescribed generalized modal-mass ratio.', ...
    'shapeNormalization','Physical displacement columns have Euclidean norm one; largest component real positive.');
if ~valid,return;end
[frequency,order]=sort(abs(allPoles(selected))/(2*pi)); selected=selected(order);
lambda=allPoles(selected); Phi=V(1:2,selected); reference=zeros(1,2);
for j=1:2
    Phi(:,j)=Phi(:,j)/norm(Phi(:,j));
    [~,reference(j)]=max(abs(Phi(:,j)));
    Phi(:,j)=Phi(:,j)*exp(-1i*angle(Phi(reference(j),j)));
end
model.lambda=lambda;
model.theta=[frequency.',(-real(lambda)./abs(lambda)).'];
model.displacementShapes=Phi;
model.shapeReference=reference;
model.shapeRatios=Phi./Phi(1,:);
model.modeEquationRelativeError= ...
    norm(M*Phi*diag(lambda.^2)+C*Phi*diag(lambda)+K*Phi,'fro') / ...
    ((norm(M,2)*max(abs(lambda))^2+norm(C,2)*max(abs(lambda))+norm(K,2))*norm(Phi,'fro'));
end
