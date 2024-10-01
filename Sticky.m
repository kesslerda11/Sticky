% This simulates Nsp species, with imigration and env. noise
% NO demographic stoch
Nsp=200;
K0=10000;
K=Nsp*K0;
pop=K0*ones(Nsp,1);
r0=1;
r_birth=r0*ones(Nsp,1);
%t_env=10/r0;
t_env=8;
%sig_env=0.025;
sig_env=0.025;
rb_min=max(0,r_birth*(1-3*sig_env*sqrt(t_env)));
rb_max=2*r0-rb_min;
t_ext=cell(Nsp,1);
for isp=1:Nsp
    t_ext{isp}=[];
end
dt=0.005;
t=0;
tf=800000;
totpop=sum(pop);
%r_imm=3;
r_imm=0.6;
%dtout=t_env;
dtout=floor(t_env/10/dt)*dt;
%tout=t_env;
tout=tf/2;
iout=0;
popt=zeros(Nsp,floor(tf/dtout/2));
zipf=zeros(Nsp,1);
lzipf=zeros(Nsp,1);
sphft=zeros(floor(tf/dtout/2),1);
rb1t=zeros(floor(tf/dtout/2),1);
while (t<tf)
    t=t+dt;
    % Assign new birth rates
    r_birth=max(rb_min,r0+(r_birth-r0)*exp(-dt/t_env)+sig_env*sqrt(2*dt)*randn(Nsp,1));
    r_birth=min(r_birth,rb_max);
    % Births
    %births=poidev(dt*r_birth.*pop);
    births=dt*r_birth.*pop;
    % Deaths
    %deaths=bnldev(pop,r0*dt*totpop/K*ones(Nsp,1));
    deaths=dt*r0*totpop/K*pop;
    % Immigratio
    %imms=poidev(dt*r_imm*ones(Nsp,1));
    imms=dt*r_imm*ones(Nsp,1);
    pop=pop + births - deaths + imms;
    if (min(pop)<0)
        disp(min(pop))
        pop=max(0,pop);
    end
    totpop=sum(pop);
    % Store extinction times
    %x_ndxs=find(pop==0);
    % for i=1:length(x_ndxs)
    %     ndx=x_ndxs(i);
    %     if isempty(t_ext{ndx})
    %         t_ext{ndx}=t;
    %     elseif (t-t_ext{ndx}(end)>t_env)
    %         t_ext{ndx}=[t_ext{ndx},t];
    %     end
    % end
    if (t>tout)
        % spop=sort(pop);
        % semilogy(max(spop,.1),'o');
        % zipf=zipf+spop;
        % lzipf=lzipf+log(max(spop,0.1));
        % title('t: ',num2str(t));
        tout=tout + dtout;
        iout=iout + 1;
        popt(:,iout)=pop;
        % drawnow;
        % Calculate Shalf
        spop=sort(pop,'descend');
        cs=cumsum(spop);
        ndx=find(cs>totpop/2,1,'first');
        if (ndx==1)
            frac=(totpop/2)/(cs(ndx));
        else
            frac=(totpop/2-cs(ndx-1))/(cs(ndx)-cs(ndx-1));
        end
        sphft(iout)=ndx-1+frac;
        rb1t(iout)=r_birth(1);
        25;
    end
end
mean(sum(popt(:,floor(end/3):end),1))/K-1
mean(sum(popt(:,floor(end/2):end),1))/K-1
mean(sum(popt(:,floor(2*end/3):end),1))/K-1
% zipf=zipf/iout;
% lzipf=lzipf/iout;
% maxpop=max(popt(1:end));
% duration=cell(ceil(log2(maxpop)),1);
% for isp=1:Nsp
%     for i=2:length(t_ext{isp})
%         ti=t_ext{isp}(i);
%         ndx=floor(t_ext{isp}(i)/dtout);
%         while (ndx>ceil(t_ext{isp}(i-1)/dtout))
%             dndx=max(1,round(log2(popt(isp,ndx))));
%             duration{dndx}=[duration{dndx},ti-ndx*dtout];
%             ndx=ndx-1;
%         end
%     end
% end
rd1t=sum(popt,1)'/K-1;
vard=var(rd1t)
net=rb1t-rd1t;
net=net-mean(net);
[xcb,lags]=xcorr(rb1t-mean(rb1t),50);
[xcd,lags]=xcorr(rd1t-mean(rd1t),50);
[xcn,lags]=xcorr(net,50);
Sb=var(rb1t)*trapz(lags*dtout,xcb/xcb(51))
Sd=var(rd1t)*trapz(lags*dtout,xcd/xcd(51))
Snet=var(net)*trapz(lags*dtout,xcn/xcn(51))
[lhy,lhx]=hist(log(popt(1:end)),100);
lpy=log(exp(-lhx).*lhy/trapz(exp(lhx),exp(-lhx).*lhy));
