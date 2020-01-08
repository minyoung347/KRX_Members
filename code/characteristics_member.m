
%% Calculating Trending 
% Corr of ratio & net_buy with price_return

member_corr_net = cell(10,75);
member_corr_ratio = cell(10,75);

for f = 1:10
    sprintf('%d',f)
    for i = (f-1)*121+1 : f*121
        for y = 1:11
            for ii = 1:75
                if sum(sign(inv_mem_net{i,2}(year_day{y,1},ii))) > size(year_day{y,1},1)/1.5
                    % for net buy
                    member_corr_net{f,ii} = [member_corr_net{f,ii}; corr(price_return{i,1}(year_day{y,1}), inv_mem_net{i,1}(year_day{y,1},ii))];

                    % for ratio
                    temp = inv_mem_net{i,1}(year_day{y,1},ii)./inv_mem_net{i,2}(year_day{y,1},ii);
%                     [res, ind] = delniz(temp,[1,2],1);
%                     if length(res) > 3
%                         member_corr_ratio{f,ii} = [member_corr_ratio{f,ii}; corr(price_return{i,1}(year_day{y,1}(ind)), res)];
%                     end
                    temp(isnan(temp)) = 0;
                    temp(isinf(temp)) = 0;
                    member_corr_ratio{f,ii} = [member_corr_ratio{f,ii}; corr(price_return{i,1}(year_day{y,1}), temp)];
                end
            end
        end
    end
end

%% Correlation between inventory according to trending

% mean & length over 121 firms and 11 years
for i = 1:41
    domestic_corr_mean(i,1) = mean(member_corr_ratio{1,domestic(i)});
    domestic_corr_length(i,1) = length(member_corr_ratio{1,domestic(i)});
end
for i = 1:21
    foreign_corr_mean(i,1) = mean(member_corr_ratio{1,foreign(i)});
    foreign_corr_length(i,1) = length(member_corr_ratio{1,foreign(i)});
end

domestic_remain = domestic(domestic_corr_length >= 100);
foreign_remain = foreign(foreign_corr_length >= 100);


[df_sort, df_sort_i] = sort([domestic_corr_mean;foreign_corr_mean]);
df = [domestic,foreign];
df_sort_ind = df(df_sort_i);

df_sort_remain = [];
for i = 1:size(df_sort_ind,2)
    if sum(domestic_remain == df_sort_ind(i)) || sum(foreign_remain == df_sort_ind(i))
        df_sort_remain = [df_sort_remain; df_sort_ind(i)];
    end
end

% 'a': the averaged correlation over 1 decile and 11 years
a = zeros(75,75);
for i = 1:121
    for y = 1:11
        temp = corr(inv_mem_net{i,1}(year_day{y,1},:));
        temp(isnan(temp)) = 0;
        temp(isinf(temp)) = 0;
        a = a + temp;
    end
end

for i = 1:75
    a(i,i) = 0;
end

figure
set(gcf,'color','w')
imagesc(a(df_sort_remain,df_sort_remain)/(121*11))
colormap(jet)
set(gca,'xtick',1:length(df_sort_remain),'xticklabel',inv_mem_name(df_sort_remain),'fontsize',9)
xtickangle(90)
set(gca,'ytick',1:length(df_sort_remain),'yticklabel',inv_mem_name(df_sort_remain),'fontsize',9)
colorbar
title('Domestic + Foreign')





[d_sort, d_sort_i] = sort(domestic_corr_mean);
d_sort_ind = domestic(d_sort_i);

d_sort_remain = [];
for i = 1:size(d_sort_ind,2)
    if sum(domestic_remain == d_sort_ind(i))
        d_sort_remain = [d_sort_remain; d_sort_ind(i)];
    end
end

[f_sort, f_sort_i] = sort(foreign_corr_mean);
f_sort_ind = foreign(f_sort_i);

f_sort_remain = [];
for i = 1:size(f_sort_ind,2)
    if sum(foreign_remain == f_sort_ind(i))
        f_sort_remain = [f_sort_remain; f_sort_ind(i)];
    end
end



figure
set(gcf,'color','w')

subplot(1,2,1)
imagesc(a(d_sort_remain,d_sort_remain)/(121*11))
colormap(jet)
set(gca,'xtick',1:length(d_sort_remain),'xticklabel',inv_mem_name(d_sort_remain),'fontsize',9)
xtickangle(90)
set(gca,'ytick',1:length(d_sort_remain),'yticklabel',inv_mem_name(d_sort_remain),'fontsize',9)
caxis([-.14 .2])
colorbar
title('Domestic')

subplot(1,2,2)
imagesc(a(f_sort_remain,f_sort_remain)/(121*11))
colormap(jet)
set(gca,'xtick',1:length(f_sort_remain),'xticklabel',inv_mem_name(f_sort_remain),'fontsize',9)
xtickangle(90)
set(gca,'ytick',1:length(f_sort_remain),'yticklabel',inv_mem_name(f_sort_remain),'fontsize',9)
caxis([-.14 .2])
colorbar
title('Foreign')

%% calculating directionality

ratio_dist = cell(10,75);
for f = 1:10
    for i = (f-1)*121+1 : f*121
        for ii = 1:75
            temp = firm_net{ii,1}(:,i)./firm_net{ii,2}(:,i);
            temp = delniz(temp,[1,2],0);
            ratio_dist{f,ii} = [ratio_dist{f,ii}; temp];
        end
    end
end



threshold = .1:.1:.9;
% lr_indv = zeros(10,length(threshold));
% lr_fore = zeros(10,length(threshold));
% lr_inst = zeros(10,length(threshold));
directionality_member = cell(75,1);
for i = 1:75
    directionality_member{i,1} = zeros(10,length(threshold));
end
for f = 1:10
    for th = 1:length(threshold)
%         for ii = 1:3
%             if ii == 1
%                 lr_indv(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
%             elseif ii == 2
%                 lr_fore(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
%             elseif ii == 3
%                 lr_inst(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
%             end
%         end
        for ii = 1:75
            directionality_member{ii,1}(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
        end
    end
end
clear ratio_dist


figure
set(gcf,'color','w')
subplot(1,3,1)
plot(directionality_member{1,1})
ylim([0 1])
xlabel('percentile')
ylabel('portion of ratio > th')
title('indv')
subplot(1,3,2)
plot(directionality_member{2,1})
ylim([0 1])
xlabel('percentile')
ylabel('portion of ratio > th')
title('fore')
subplot(1,3,3)
plot(directionality_member{3,1})
ylim([0 1])
xlabel('percentile')
ylabel('portion of ratio > th')
title('inst')


%% directionality & trending for members

for i = 1:10
    for ii = 1:75
        trending(i,ii) = mean(delniz(member_corr_ratio{i,ii},[1,2],0));
%         trending(i,ii) = mean(delniz(member_corr_net{i,ii},[1,2],0));
    end
end

k = 1;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
color_3 = linspecer(3);



% directionality & trending for 41 domestic and 21 foreign
figure2 = figure;
set(gcf,'color','w')
subplot(1,2,1)
for ii = 1:41
    plot(trending(:,domestic(ii)),directionality_member{domestic(ii),1}(:,k),'--k','linewidth',1.5)
    hold on
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
% xlim([-.5 .5])
% ylim([.4 1])
for ii = 1:41
    for i = 1:10
        scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
        hold on
    end
end

subplot(1,2,2)
for ii = 1:21
    plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'b','linewidth',1.5)
    hold on
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
for ii = 1:21
    for i = 1:10
        scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
%         scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',[0 0 1])
        hold on
    end
end
% text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.6 .35])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
lgd1 = legend('Individual','Foreigner','Institution');
lgd1.FontSize = 12;
legend boxoff



%% figure in the paper

figure1 = figure;
set(gcf,'color','w')

subplot(1,3,1)
plot(trending(:,1),directionality_member{1,1}(:,k),'-','color',color_3(:,3),'linewidth',1.5)
hold on
plot(trending(:,2),directionality_member{2,1}(:,k),'--','color',color_3(:,2),'linewidth',1.5)
hold on
plot(trending(:,3),directionality_member{3,1}(:,k),':','color',color_3(:,1),'linewidth',1.5)
hold on
plot([0 0],[0 1],'k--')
hold on
plot([-.6 .35],[.5 .5],'k--')
for i = 1:10
    hold on
    scatter(trending(i,1),directionality_member{1,1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
    hold on
    scatter(trending(i,2),directionality_member{2,1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
    hold on
    scatter(trending(i,3),directionality_member{3,1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
end 
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.6 .35])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
lgd1 = legend('Individual','Foreigner','Institution');
lgd1.FontSize = 12;
legend boxoff

subplot(1,3,2)
for ii = 1:41
    if ca(ii,3) < -.00 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 1:1
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(3,:))
            hold on
        end
        text(trending(1,domestic(ii))+.01,directionality_member{domestic(ii),1}(1,k),inv_mem_name(domestic(ii)),'fontsize',7)
    elseif ca(ii,3) > .00 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 1:1
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(1,:))
            hold on
        end
        text(trending(1,domestic(ii))+.01,directionality_member{domestic(ii),1}(1,k),inv_mem_name(domestic(ii)),'fontsize',7)
    end
end

for ii = 1:21
    if cb(ii,2) > .1 && mean_vol_member(foreign(ii)) > mean(mean_vol_member(foreign))/10
%         plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 1:1
            scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(2,:))
            hold on
        end
        text(trending(1,foreign(ii))+.01,directionality_member{foreign(ii),1}(1,k),inv_mem_name(foreign(ii)),'fontsize',7)
    end
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.45 .28])
ylim([.43 1])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
box on



subplot(1,3,3)
for ii = 1:41
    if ca(ii,3) < -.02 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 10:10
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(3,:))
            hold on
        end
        text(trending(10,domestic(ii))+.01,directionality_member{domestic(ii),1}(10,k),inv_mem_name(domestic(ii)),'fontsize',7)
    elseif ca(ii,3) > .02 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 10:10
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(1,:))
            hold on
        end
        text(trending(10,domestic(ii))+.01,directionality_member{domestic(ii),1}(10,k),inv_mem_name(domestic(ii)),'fontsize',7)
    end
end

for ii = 1:21
    if cb(ii,2) > .1 && mean_vol_member(foreign(ii)) > mean(mean_vol_member(foreign))/10
%         plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 10:10
            scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(2,:))
            hold on
        end
        text(trending(10,foreign(ii))+.01,directionality_member{foreign(ii),1}(10,k),inv_mem_name(foreign(ii)),'fontsize',7)
    end
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
text(-0.13,1.08,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.45 .28])
ylim([.43 1])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
box on



dx=0.06;
dy=0.06;
x=(1-4*dx)/2.9;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(3),'Position',[2*dx-dx,3*dy,x,y]);
set(AxesHandle(2),'Position',[x+2.6*dx-dx/2,3*dy,x,y]);
set(AxesHandle(1),'Position',[2*x+3.5*dx-dx/3,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1800, 600])







