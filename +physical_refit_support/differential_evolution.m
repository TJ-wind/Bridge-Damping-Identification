function [q,value,flag,output,points,history]=differential_evolution(fun,population,maxIterations,maxStallIterations)
%DIFFERENTIAL_EVOLUTION Original bounded rand/1/bin fallback, F=.8, CR=.9.
[n,d]=size(population); values=zeros(n,1);
for j=1:n,values(j)=fun(population(j,:));end
[value,idx]=min(values); q=population(idx,:);
history=[0,n,value,mean(values),q]; count=n; stall=0; flag=0;
for iteration=1:maxIterations
    old=value;
    for j=1:n
        pool=[1:j-1,j+1:n]; draw=pool(randperm(n-1,3));
        donor=population(draw(1),:)+.8*(population(draw(2),:)-population(draw(3),:));
        donor=min(max(donor,0),1); mask=rand(1,d)<.9; mask(randi(d))=true;
        trial=population(j,:); trial(mask)=donor(mask);
        candidate=fun(trial); count=count+1;
        if candidate<values(j),population(j,:)=trial;values(j)=candidate;end
    end
    [value,idx]=min(values); q=population(idx,:);
    history(end+1,:)=[iteration,count,value,mean(values),q]; %#ok<AGROW>
    if old-value<=1e-8*max(1,abs(old)),stall=stall+1;else,stall=0;end
    if stall>=maxStallIterations,flag=1;break;end
end
output=struct('iterations',iteration,'funccount',count, ...
    'message',sprintf('Bounded DE rand/1/bin; F=.8, CR=.9; at most %d generations, %d-stall rule.',maxIterations,maxStallIterations));
points=struct('X',population,'Fval',values);
end
