
%% correlation between inventory variation & return (Distributions)
% delay에 따라서 의미 달라짐 (delay = -1, 0, 1)

% cir: Correlation b/w Iventory & Return
cir_1 = cell(10,75);
cir_2 = cell(10,75);
cir_3 = cell(10,75);

for ii = 1:75
    sprintf('%d',ii)
    for f = 1:10
        temp_1 = [];
        temp_2 = [];
        temp_3 = [];
        for i = 1+121*(f-1):121*f
            for y = 1:11
                if length(inv_mem_net{i,1}(year_day{y,1},ii)) > length(year_day{y,1})/1.5
                    temp_1 = [temp_1; corr(inv_mem_net{i,1}(year_day{y,1},ii), price_return{i,1}(year_day{y,1},1))];
                    temp_2 = [temp_2; corr(inv_mem_net{i,1}(year_day{y,1}(1:end-1),ii), price_return{i,1}(year_day{y,1}(2:end),1))];
                    temp_3 = [temp_3; corr(inv_mem_net{i,1}(year_day{y,1}(2:end),ii), price_return{i,1}(year_day{y,1}(1:end-1),1))];
                end
            end
        end
        cir_1{f,ii} = temp_1;
        cir_2{f,ii} = temp_2;
        cir_3{f,ii} = temp_3;
    end
end

% reversed: cir_1{f,ii}, cc_1(ii,f)
for ii = 1:75
    for f = 1:10
        cc_1(ii,f) = mean(delniz(cir_1{f,ii},[1,2],0));
        cc_2(ii,f) = mean(delniz(cir_2{f,ii},[1,2],0));
        cc_3(ii,f) = mean(delniz(cir_3{f,ii},[1,2],0));
    end
end


figure
set(gcf,'color','w')
for i = 1:41
    subplot(6,7,i)
%     subplot(1,3,i)
    histfit(cir_1{1,domestic_sort(i)}) % histogram of cir for 1 decile and domestic members
    hold on
    plot([0 0],[0 100],'k','linewidth',3)
    xlim([-.5 .5])
    title(sprintf('%s, %.2d',inv_mem_name{domestic_sort(i,1),1},mean(delniz(cir_1{1,domestic_sort(i,1)},[1,2],0))))
end


figure
set(gcf,'color','w')
for i = 1:21
    subplot(3,7,i)
    histfit(cir_1{1,foreign_sort(i,2)}) % histogram of cir for 1 decile and foreign members
    hold on
    plot([0 0],[0 100],'k','linewidth',3)
    xlim([-.5 .5])
    title(sprintf('%s, %.2d',inv_mem_name{foreign_sort(i,2),1},mean(delniz(cir_1{1,foreign_sort(i,2)},[1,2],0))))
end


%% Feedback & Price impact

a = zeros(50,75);
percent1 = 0;

for i = 1:121 %1210-121+1:1210
    
    temp1 = find(price_return{i,1}(1:2722) > log((100+percent1)/100));

    if ~isempty(temp1)
        for ii = 1:50
            delay = ii-25;
            temp_d = temp1(find(temp1+delay>0 & temp1+delay<=2722));
            if ~isempty(temp_d)
                a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:).*price{i,2}(temp_d+delay,:));
%                 a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:));
            end
        end
    end
end


% plot for individuals, foreigners, institutions
figure
set(gcf,'color','w')
for i = 1:3
    subplot(1,3,i)
    plot(-24:1:25,a(:,i),'k')
    hold on
    scatter(0,a(25,i),'k','filled')
    hold on
    scatter(-1,a(24,i),'k')
    hold on
    scatter(-2,a(23,i),'k')
    hold on
    scatter(1,a(26,i),'k')
    hold on
    scatter(2,a(27,i),'k')
    title(inv_mem_name{i,1},'fontsize',18)
    xlim([-25 25])
    xlabel('Time lag (day)','fontsize',16)
    ylabel('Net price (Won)','fontsize',16)
end


% plot mean & std
figure0 = figure;
set(gcf,'color','w')
subplot(1,3,1)
ind = [];
for i = 1:41
    if ca(i,1) < -.0
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k','linewidth',1)
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,1) color_3(2,1) color_3(3,1)])
alpha(.5)
hold on
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (institution dominent)')
legend boxoff


subplot(1,3,2)
ind = [];
for i = 1:41
    if ca(i,1) > .05
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k')
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,3) color_3(2,3) color_3(3,3)])
alpha(.5)
hold on
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (individual dominent)')
legend boxoff

subplot(1,3,3)
plot(-24:1:25,mean(a(:,foreign),2),'k')
hold on
scatter(-1:1,mean(a(24:26,foreign),2),'k')
hold on
ciplot(mean(a(:,foreign),2)-std(a(:,foreign)')',mean(a(:,foreign),2)+std(a(:,foreign)')',-24:1:25,[color_3(1,2) color_3(2,2) color_3(3,2)])
alpha(.5)
hold on
text(-0.13,1.08,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Foreign members')
legend boxoff

dx=0.06;
dy=0.06;
x=(1-4*dx)/2.9;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure0,'Type','axes');
set(AxesHandle(3),'Position',[2*dx-dx,3*dy,x,y]);
set(AxesHandle(2),'Position',[x+2.6*dx-dx/2,3*dy,x,y]);
set(AxesHandle(1),'Position',[2*x+3.5*dx-dx/3,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1200, 420])



% plot the mean of all members
figure
set(gcf,'color','w')
subplot(1,3,1)
for i = 1:41
    if ca(i,3) < -.05 && mean_vol_member(domestic(i)) > mean(mean_vol_member(domestic))/5
        plot(-24:1:25,a(:,domestic(i)),'k')
        hold on
    end
end
subplot(1,3,2)
for i = 1:41
    if ca(i,3) > .1 && mean_vol_member(domestic(i)) > mean(mean_vol_member(domestic))/5
        plot(-24:1:25,a(:,domestic(i)),'k')
        hold on
    end
end
subplot(1,3,3)
for i = 1:21
    plot(-24:1:25,a(:,foreign(i)),'k')
    hold on
end


%% Feedback & Price impact
% simultaneously 1 dcile & 10 decile

a = zeros(50,75);
percent1 = 0;

for i = 1:121
    
    temp1 = find(price_return{i,1}(1:2722) > log((100+percent1)/100));

    if ~isempty(temp1)
        for ii = 1:50
            delay = ii-25;
            temp_d = temp1(find(temp1+delay>0 & temp1+delay<=2722));
            if ~isempty(temp_d)
                a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:).*price{i,2}(temp_d+delay,:));
%                 a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:));
            end
        end
    end
end





figure0 = figure;
set(gcf,'color','w')


subplot(2,3,1)
ind = [];
for i = 1:41
    if ca(i,1) < -.0
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k','linewidth',1)
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,1) color_3(2,1) color_3(3,1)])
alpha(.5)
hold on
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (institution dominent)')
legend boxoff


subplot(2,3,2)
ind = [];
for i = 1:41
    if ca(i,1) > .05
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k')
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,3) color_3(2,3) color_3(3,3)])
alpha(.5)
hold on
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (individual dominent)')
legend boxoff

subplot(2,3,3)
plot(-24:1:25,mean(a(:,foreign),2),'k')
hold on
scatter(-1:1,mean(a(24:26,foreign),2),'k')
hold on
ciplot(mean(a(:,foreign),2)-std(a(:,foreign)')',mean(a(:,foreign),2)+std(a(:,foreign)')',-24:1:25,[color_3(1,2) color_3(2,2) color_3(3,2)])
alpha(.5)
hold on
text(-0.13,1.08,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Foreign members')
legend boxoff



a = zeros(50,75);
percent1 = 0;

for i = 1210-121+1:1210
    
    temp1 = find(price_return{i,1}(1:2722) > log((100+percent1)/100));

    if ~isempty(temp1)
        for ii = 1:50
            delay = ii-25;
            temp_d = temp1(find(temp1+delay>0 & temp1+delay<=2722));
            if ~isempty(temp_d)
                a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:).*price{i,2}(temp_d+delay,:));
%                 a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:));
            end
        end
    end
end


subplot(2,3,4)
ind = [];
for i = 1:41
    if ca(i,1) < -.0
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k','linewidth',1)
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,1) color_3(2,1) color_3(3,1)])
alpha(.5)
hold on
text(-0.13,1.08,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (institution dominent)')
legend boxoff


subplot(2,3,5)
ind = [];
for i = 1:41
    if ca(i,1) > .05
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k')
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,3) color_3(2,3) color_3(3,3)])
alpha(.5)
hold on
text(-0.13,1.08,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (individual dominent)')
legend boxoff

subplot(2,3,6)
plot(-24:1:25,mean(a(:,foreign),2),'k')
hold on
scatter(-1:1,mean(a(24:26,foreign),2),'k')
hold on
ciplot(mean(a(:,foreign),2)-std(a(:,foreign)')',mean(a(:,foreign),2)+std(a(:,foreign)')',-24:1:25,[color_3(1,2) color_3(2,2) color_3(3,2)])
alpha(.5)
hold on
text(-0.13,1.08,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Foreign members')
legend boxoff


dx=0.05;
dy=0.05;
x=(1-5*dx)/2.9;
y=(1-5*dy)/1.9;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure0,'Type','axes');
set(AxesHandle(6),'Position',[dx,y+3.45*dy,x,y]);
set(AxesHandle(5),'Position',[x+2.5*dx,y+3.45*dy,x,y]);
set(AxesHandle(4),'Position',[2*x+4*dx,y+3.45*dy,x,y]);
set(AxesHandle(3),'Position',[dx,dy+dy/2,x,y]);
set(AxesHandle(2),'Position',[x+2.5*dx,dy+dy/2,x,y]);
set(AxesHandle(1),'Position',[2*x+4*dx,dy+dy/2,x,y]);
set(gcf,'Position',[100, 100, 1300, 750])

