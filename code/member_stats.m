
%% member trading on/off
% output: [start_end, start_end_ind]
% start_end: start_day, end_day, # of days
% start_end_ind: index of start_end

total_volume = cell(1210,1);
agg_total_volume = zeros(2722,75);
for f = 1:1210
    f
    total_volume{f,1} = zeros(2722,75);
    for i = 1:2722
        for ii = 1:75
            total_volume{f,1}(i,ii) = sign(firm_net{ii,2}(i,f));
        end
    end
    
    agg_total_volume = sign(agg_total_volume + total_volume{f,1});
end


f = 1;
figure
set(gcf,'color','w')
subplot(1,2,1)
imagesc(agg_total_volume(:,[domestic,foreign]))
set(gca,'xtick',1:62,'xticklabel',inv_mem_name([domestic,foreign]),'fontsize',10)
xtickangle(90)
title('All firm')
subplot(1,2,2)
imagesc(total_volume{f,1}(:,[domestic,foreign]))
set(gca,'xtick',1:62,'xticklabel',inv_mem_name([domestic,foreign]),'fontsize',10)
xtickangle(90)
title(sprintf('%s',name_firm{f,1}))

figure
set(gcf,'color','w')
for i = 1:40
    subplot(4,10,i)
    imagesc(total_volume{i+0,1}(:,[domestic,foreign]))
    title(sprintf('%s',name_firm{i+0,1}))
end

start_end = zeros(75,3);
start_end_ind = zeros(75,2);
for i = 1:75
    start_end(i,1) = day{1,1}(min(find(agg_total_volume(:,i) == 1)));
    start_end(i,2) = day{1,1}(max(find(agg_total_volume(:,i) == 1)));
    start_end(i,3) = length(find(agg_total_volume(:,i) == 1));
    
    start_end_ind(i,1) = min(find(agg_total_volume(:,i) == 1));
    start_end_ind(i,2) = max(find(agg_total_volume(:,i) == 1));
end

clear total_volume agg_total_volume

%% statistics for member firms
% Market capitalization 1 decile stocks에 대해서, member trading 통계
% output: [total_pvol, active_day, active_stock]
% total_pvol: 전체기간에 대해서 회원사가 1 decile stocks을 거래한 pvolume
% active_day: 회원사가 활동한 기간
% active_stock: 회원사가 거래한 주식 수

total_pvol = zeros(1,75);
active_day = zeros(1,75);
active_stock = zeros(1,75);

for ii = 1:75
    active_temp = zeros(2722,121);
    for i = 1 : 121
        pvol = inv_mem_net{i,2}(1:2722,ii) .* price{i,2}(1:2722);
        total_pvol(1,ii) = total_pvol(1,ii) + sum(pvol);
        
        active_temp(:,i) = inv_mem_net{i,2}(1:2722,ii)>0;
    end
    active_stock(1,ii) = length(find(sum(active_temp,1) > 0));
    active_day(1,ii) = length(find(sum(active_temp,2) > 0));
end

figure
set(gcf,'color','w')
subplot(1,2,1)
bar(total_pvol(domestic))
ylim([0 8*10^15])
set(gca,'xtick',1:41,'xticklabel',inv_mem_name(domestic),'fontsize',9)
xtickangle(90)
subplot(1,2,2)
bar(total_pvol(foreign))
ylim([0 8*10^15])
set(gca,'xtick',1:21,'xticklabel',inv_mem_name(foreign),'fontsize',9)
xtickangle(90)

figure
set(gcf,'color','w')
subplot(1,2,1)
bar(total_pvol(domestic)./active_day(domestic))
ylim([0 3*10^12])
set(gca,'xtick',1:41,'xticklabel',inv_mem_name(domestic),'fontsize',9)
xtickangle(90)
subplot(1,2,2)
bar(total_pvol(foreign)./active_day(foreign))
ylim([0 3*10^12])
set(gca,'xtick',1:21,'xticklabel',inv_mem_name(foreign),'fontsize',9)
xtickangle(90)

figure
set(gcf,'color','w')
subplot(1,2,1)
bar(total_pvol(domestic)./active_day(domestic)./active_stock(domestic))
ylim([0 2.5*10^9])
set(gca,'xtick',1:41,'xticklabel',inv_mem_name(domestic),'fontsize',9)
xtickangle(90)
subplot(1,2,2)
bar(total_pvol(foreign)./active_day(foreign)./active_stock(foreign))
ylim([0 2.5*10^9])
set(gca,'xtick',1:21,'xticklabel',inv_mem_name(foreign),'fontsize',9)
xtickangle(90)

a2 = inv_mem_name([domestic,foreign]);
a3 = start_end([domestic,foreign],1:2);
a3p = active_day([domestic,foreign]);
a4 = active_stock([domestic,foreign]);
a5 = total_pvol([domestic,foreign]);
a6 = total_pvol([domestic,foreign])./active_day([domestic,foreign])./active_stock([domestic,foreign]);

fid = fopen('stat_member.csv','w');
for i = 1:62
    start_ym = num2str(a3(i,1));
    end_ym = num2str(a3(i,2));
    fprintf(fid,'%.0f,%s,%s-%s (%.0f),%.0f,%.2f,%.0f\n',i,a2{i,1},start_ym(3:6),end_ym(3:6),a3p(i),a4(i),a5(i)/10^12,a6(i)/10^6);
end
fclose(fid);


