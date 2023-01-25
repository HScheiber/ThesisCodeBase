epsilon = 1;
sigma = 1;
%gammas_cell = {[0.1:0.01:0.3 0.35:0.05:5 5.1:0.002:5.9] [6.9:-0.05:6.45 6.4:-0.002:6.1] 7.0 7.0:0.01:8.5};
gammas = 3.3:0.01:20;

title_txt = {'$0 \leq \gamma \leq 20$'};
rmin = 0;
rmax = 3;
ylims = [-2 5];
fs=32;

% Update plot
figh = figure('WindowState','Maximized');
T = tiledlayout(figh,1,1);

%yline(ax,0,':k',"LineWidth",1)


ax = nexttile(T);
hold(ax,'on')

Colours = cbrewer('seq','Blues',max([numel(gammas) 3]),'spline');

for jdx = 1:numel(gammas)
    gamma = gammas(jdx);

    ffun = @(r) (2*epsilon/(1 - (3/(gamma+3)))).*((sigma^6)).*((3/(gamma+3)).*exp(gamma.*(1 - (r./sigma))) - 1./((sigma^6) + (r.^6)));
    fpfun = @(r) -(2*epsilon/(1 - (3/(gamma+3)))).*((6*sigma^6).*(r.^5))./((sigma^6 + (r.^6)).^2).*((3/(gamma+3)).*exp(gamma.*(1 - (r./sigma))) - 1) ...
        - (2*epsilon/(1 - (3/(gamma+3)))).*((sigma^6)./((sigma^6) + (r.^6))).*(3/(gamma+3)).*((gamma.*exp(gamma.*(1 - r./sigma)))./sigma);
    
    rt = fzero(fpfun,sigma);
    fmin = ffun(rt);
    
    r = rmin:0.001:rmax;
    xmax = rmax;

    % Calculate potential
    f = ffun(r);

    % Derivative of potential
    fprime = fpfun(r);

    if numel(gammas) == 1
        plot(ax,r,f,'-b',"LineWidth",3,'color',Colours(end,:))
    else
        plot(ax,r,f,'-b',"LineWidth",3,'color',Colours(jdx,:))
    end
end
plot(ax,[0 xmax],[0 0],'--k',"LineWidth",1)

plot(ax,[sigma sigma],[ffun(sigma) ylims(1)],':r',"LineWidth",3)
plot(ax,[rt rt],[ffun(rt) ylims(1)],':r',"LineWidth",3)

plot(ax,[0 sigma],[ffun(sigma) ffun(sigma)],':r',"LineWidth",3)
plot(ax,[0 rt],[fmin fmin],':r',"LineWidth",3)

scatter(ax,rt,fmin,80,'r','filled','o',"LineWidth",1,'MarkerEdgeColor','k')
scatter(ax,sigma,ffun(sigma),80,'r','filled','o',"LineWidth",1,'MarkerEdgeColor','k')
ylim(ax,ylims)
xlim(ax,[rmin,xmax])
title(ax,title_txt,'fontsize',fs,'Interpreter','latex')

set(ax,'box','on','TickLabelInterpreter','latex');
set(ax,'XMinorTick','off','YMinorTick','off','FontSize',fs);
xlabel(ax,'$r$','fontsize',fs,'Interpreter','latex');

%ylabel(ax,'$u^{BH}$','fontsize',fs,'Interpreter','latex');

[xtt,xidx] = sort([0 sigma]);
xlab = {'0' '$\sigma$'};
xlabsrt = xlab(xidx);

[ytt,yidx] = sort([-epsilon 0]);
ylab = {'$-\varepsilon$' '0'};
ylabsrt = ylab(yidx);

xticks(ax,xtt)
yticks(ax,ytt)
xticklabels(ax,xlabsrt)
yticklabels(ax,ylabsrt)

%exportgraphics(figh,'WBK_pot.eps')
%exportgraphics(figh,'WBK_pot.png','resolution',600)