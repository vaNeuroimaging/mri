function prettyscatterplot(xvals,yvals,fitline)

figure('Color','white','position',[1982 478 1352 804]);
plot(xvals,yvals,'k.','MarkerSize',40)
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
box off

if fitline
    hold on
    FO = fit(xvals(:),yvals(:),'poly1');
    
    extendamount = range(xvals) ./ 20;
    
    fittedxvals = [(min(xvals) -extendamount) ; sort(xvals(:),'ascend') ; (max(xvals)+extendamount)];
    fittedyvals = (FO.p1 .* fittedxvals) + FO.p2;
    
    plot(fittedxvals,fittedyvals,'k-','LineWidth',3)
    xlim([(min(xvals) -extendamount) (max(xvals) +extendamount)]);
end