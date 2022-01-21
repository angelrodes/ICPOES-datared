function [ out ] = linear_regression_chisq_fn(x,dx,y,dy)
% Perform a linear fit to data with x and y uncertainties (one-sigma)
% Outputs: Output arrays are oredered from best to worst fitting within 1-sigma
%   slopes: array with possible sloples within 1 sigma confidence
%   intercepts: array with possible intrcepts within 1 sigma confidence
%     One-sigma fitting lines will be y = x .* slopes + intercepts
%   chisq: array of chi-squared values corresponding to the slopes and itercepts values
%   redchisq: reduced chi-squared value of the best fit = minimum chisq value / degrees of freedom
%   R2: R2 of the best fit
%   sloping: 1 if slope is not 0 within the one-sigma range
% Angel Rodes, 2019

%% select number of models
testpointstotal=1e6;
testpoints=round(testpointstotal.^0.5);

%% doplot=1 if you want to see plots (for testing purposes)
doplot=0;

%% All inputs to single rows
x=x(:)';
dx=dx(:)';
y=y(:)';
dy=dy(:)';

%% Check data
check=std([length(x),length(y),length(dx),length(dy)])>0 | length(x)<3 | min(dy)<0 | min(dx)<0;
if check
    error('Check your data!')
end

%% Fix data if there are points without any uncertainty
dytemp=dy;
dxtemp=dx;
if sum(dy==0 & dx==0)>0
    % Consider minimum uncertainties or SDOM
    if sum(dy>0)==0
        dytemp(dy==0 & dx==0)=std(y)/length(y);
    else
        dytemp(dy==0 & dx==0)=min(dy(dy>0));
    end
    if sum(dx>0)==0
        dxtemp(dy==0 & dx==0)=std(x)/length(x);
    else
        dxtemp(dy==0 & dx==0)=min(dx(dx>0));
    end
end
dy=dytemp;
dx=dxtemp;

%% first random guess
st=zeros(1,testpoints);
yit=zeros(1,testpoints);
for n=1:testpoints
    xtemp=normrnd_BoxMuller(x,dx);
    ytemp=normrnd_BoxMuller(y,dy);
    slopetemp=(ytemp-mean(ytemp))/(xtemp-mean(xtemp));
    intercepttemp=mean(ytemp-slopetemp*xtemp);
    st(n)=slopetemp;
    yit(n)=intercepttemp;
end



iteration=0;
convergence=0;
while convergence<1
    iteration=iteration+1;
    %% burld matrices of slopes and intercepts
    if iteration==1
        expandfactor=10;
        sm=linspace(mean(st)-std(st)*expandfactor,mean(st)+std(st)*expandfactor,testpoints);
        yim=linspace(mean(yit)-std(yit)*expandfactor,mean(yit)+std(yit)*expandfactor,testpoints);
        [SM,YIM] = meshgrid(sm,yim);
    else
        expandfactor=max(1.5,10/iteration);
        sm=linspace(SM(selbestchi)-max(std(st),max(abs(SM(selbestchi)-SM(selchi))))*expandfactor,...
            SM(selbestchi)+max(std(st),max(abs(SM(selbestchi)-SM(selchi))))*expandfactor,...
            testpoints);
        yim=linspace(YIM(selbestchi)-max(std(yit),max(abs(YIM(selbestchi)-YIM(selchi))))*expandfactor,...
            YIM(selbestchi)+max(std(yit),max(abs(YIM(selbestchi)-YIM(selchi))))*expandfactor,...
            testpoints);
        [SM,YIM] = meshgrid(sm,yim);
    end
    %% calculate Chi-square matrix
    CHISQ=0*SM;
    for n=1:length(x)
        DX=abs( (x(n)-(y(n)-YIM)./SM) );
        DY=abs( (y(n)-(YIM+x(n)*SM)) );
        Distance=DY./(dy(n)^2+(DY*dx(n)./DX).^2).^0.5;
        CHISQ=CHISQ+Distance.^2;
    end
    
    %% calculate one sigma limits limits and goodness of fit
    minchi=min(CHISQ(:));
    dof=length(x)-2;
%     maxchi=chi2inv(0.6827+chi2cdf(minchi,dof)*(1-0.6827),dof);
    maxchi=minchi+max(1,dof); % simplified to avoid using the statistic package in octave
    if isinf(maxchi)
        maxchi=minchi+dof;
    end
    selchi=(CHISQ<maxchi);
    selbestchi=find(CHISQ==minchi,1,'first');
    redchisq=minchi/dof;
    nsolutions=sum(selchi(:));
    ybf=YIM(selbestchi)+x*SM(selbestchi);
    Rfit = corrcoef(y,ybf);
    Rsqfit = Rfit(1,2).^2;
    
    %% check convergence
    if iteration>10
        convergence=1;
    end
    if nsolutions>100 && (max(selchi(1,:))|max(selchi(end,:))|max(selchi(:,1))|max(selchi(:,end)))==0
        convergence=convergence+0.5;
    end
    
end

%% Ouputs
valuestoexport=find(selchi);
chisqvalues=CHISQ(selchi);
[~,b]=sort(chisqvalues);
valuestoexport=valuestoexport(b); % sort indexes by chi order
maxnvalues=min(1000,length(valuestoexport));
b=round(linspace(1,length(valuestoexport),maxnvalues));
valuestoexport=valuestoexport(b); % reduce the values
out.slopes=SM(valuestoexport);
out.intercepts=YIM(valuestoexport);
out.chisq=CHISQ(valuestoexport);
out.redchisq=redchisq;
out.R2=Rsqfit;
% out.prob=chi2pdf(out.chisq,dof);
out.sloping=( length(unique(sign(out.slopes)))==1 );


%% plot slope-intercept results
if doplot
    figure
    subplot(1,2,1)
    hold on
    plot(SM(:),YIM(:),'.','Color',[0.7 0.7 0.7])
    plot(SM(selchi),YIM(selchi),'.r')
    plot(SM(selbestchi),YIM(selbestchi),'.k')
    plot([mean(st)-std(st),mean(st)+std(st)],[mean(yit),mean(yit)],'-b')
    plot([mean(st),mean(st)],[mean(yit)-std(yit),mean(yit)+std(yit)],'-b')
    xlabel('slope')
    ylabel('intercept')
    title(['RedX2=' num2str(redchisq) ' ; R2=' num2str(Rsqfit)])
    % end
    
    %% Plot fitting
    % doplot=1;
    % if doplot
    %     figure
    subplot(1,2,2)
    hold on
    
    
    % plot rand fit
    plotpoints=1000;
    xplot=linspace(min(x-3*dx),max(x+3*dx),plotpoints);
    yplotbest=mean(yit)+xplot*mean(st);
    yplotmax=yplotbest;
    yplotmin=yplotbest;
    for theta=linspace(0,2*pi,plotpoints)
        yplotmax=max(yplotmax,mean(yit)+std(yit)*sin(theta)+xplot*(mean(st)+std(st)*cos(theta)));
        yplotmin=min(yplotmin,mean(yit)+std(yit)*sin(theta)+xplot*(mean(st)+std(st)*cos(theta)));
    end
    plot(xplot,yplotbest,':b')
%      plot(xplot,yplotmax,':b')
%      plot(xplot,yplotmin,':b')
    
    
    % plot samples
    for n=1:length(x)
        plotpoints=100;
        theta=linspace(0,2*pi,plotpoints);
        xplot=x(n)+dx(n)*sin(theta);
        yplot=y(n)+dy(n)*cos(theta);
        plot(xplot,yplot,'-k')
        plot(x(n),y(n),'xk')
    end
    xlabel('x')
    ylabel('y')
    
    % plot CHISQ fit
    plotpoints=1000;
    xplot=linspace(min(x-3*dx)-std(x)/length(x),max(x+3*dx)+std(x)/length(x),plotpoints);
    yplotbest=YIM(selbestchi)+xplot*SM(selbestchi);
    yplotmax=yplotbest;
    yplotmin=yplotbest;
    for n=find(selchi)'
        yplotmax=max(yplotmax,YIM(n)+xplot*SM(n));
        yplotmin=min(yplotmin,YIM(n)+xplot*SM(n));
    end
    plot(xplot,yplotbest,'-r')
    plot(xplot,yplotmax,'--r')
    plot(xplot,yplotmin,'--r')
    %     text(xplot(end),yplotbest(end),['RedX2=' num2str(redchisq) ' ; N=' num2str(nsolutions)],'Color','r')
end

end

