%% from negative herding (RMT)

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')

load('data_190314.mat');

% % add Nation data to domestic institutions
% for firm = 1:1216
%     inv_mem_net{firm,1}(1:1378,3) = inv_mem_net{firm,1}(1:1378,3) + inv_mem_net{firm,1}(1:1378,11);
%     inv_mem_net{firm,2}(1:1378,3) = inv_mem_net{firm,2}(1:1378,3) + inv_mem_net{firm,2}(1:1378,11);
% end

%% member trading on/off

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


f = 100;
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
for i = 1:75
    start_end(i,1) = day{1,1}(min(find(agg_total_volume(:,i) == 1)));
    start_end(i,2) = day{1,1}(max(find(agg_total_volume(:,i) == 1)));
    start_end(i,3) = length(find(agg_total_volume(:,i) == 1));
end


%% statistics for member firms

total_pvol = zeros(1,75);
active_day = zeros(1,75);
active_stock = zeros(1,75);

for ii = 1:75
    active_temp = zeros(2722,1210);
    for i = 1 : 1210
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



%% 
% mean trading
% sort domestic & foreign list by corr with ifd, year
% for each decile

f = 1;
num_firm = 121;
start_firm = (f-1)*num_firm + 1;
end_firm = f*num_firm;

price_ori = zeros(2722,num_firm);
for i = 1:2722
    firm_i = 1;
    for ii = start_firm:end_firm
        price_ori(i,firm_i) = price{firm_i,1}(i,1);
        firm_i = firm_i + 1;
    end
end

% mean volume of member firm
% this matrix will be used in the inventory correlation with ifd
for i = 1:size(firm_net,1)
    sum_member = sum(firm_net{i,2}(1:2722,start_firm:end_firm).*price_ori,2);
    ind = find(sum_member~=0);
    mean_vol_member(i,1) = mean(sum_member(ind));
end

figure
set(gcf,'color','w')
bar(mean_vol_member([domestic,foreign]))
set(gca,'xtick',1:62,'xticklabel',inv_mem_name([domestic,foreign]),'fontsize',9)
xlabel('member')
ylabel('mean volume')
xtickangle(90)


% 'ca', 'cb', 'cc' will be used in the inventory correlation with ifd

% 'ca' is the correlation domestic member(41) with ifd(3)
a1 = zeros(num_firm,11);
a2 = zeros(num_firm,11);
a3 = zeros(num_firm,11);
for i = 1:41
    firm_i = 1;
    for f = start_firm:end_firm
        for y = 1:11
            if length(find(inv_mem_net{f,2}(year_day{y,1},domestic(i))>0)) > length(year_day{y,1})/2
                a1(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},1),inv_mem_net{f,1}(year_day{y,1},domestic(i)));
                a2(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},2),inv_mem_net{f,1}(year_day{y,1},domestic(i)));
                a3(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},3),inv_mem_net{f,1}(year_day{y,1},domestic(i)));
            end
        end
        firm_i = firm_i + 1;
    end
    ca(i,1) = mean(delniz(a1(:),[1,3],0));
    ca(i,2) = mean(delniz(a2(:),[1,3],0));
    ca(i,3) = mean(delniz(a3(:),[1,3],0));
end


% 'cb' is the correlation foreign member(21) with ifd(3)
a1 = zeros(num_firm,11);
a2 = zeros(num_firm,11);
a3 = zeros(num_firm,11);
for i = 1:21
    firm_i = 1;
    for f = start_firm:end_firm
        for y = 1:11
            if length(find(inv_mem_net{f,2}(year_day{y,1},foreign(i))>0)) > length(year_day{y,1})/2
                a1(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},1),inv_mem_net{f,1}(year_day{y,1},foreign(i)));
                a2(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},2),inv_mem_net{f,1}(year_day{y,1},foreign(i)));
                a3(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},3),inv_mem_net{f,1}(year_day{y,1},foreign(i)));
            end
        end
        firm_i = firm_i + 1;
    end
    cb(i,1) = mean(delniz(a1(:),[1,3],0));
    cb(i,2) = mean(delniz(a2(:),[1,3],0));
    cb(i,3) = mean(delniz(a3(:),[1,3],0));
end


% 'cc' is the correlation institutions(8) with ifd(3)
a1 = zeros(num_firm,11);
a2 = zeros(num_firm,11);
a3 = zeros(num_firm,11);
ins = 4:11;
for i = 1:8
    firm_i = 1;
    for f = start_firm:end_firm
        for y = 1:11
            if length(find(inv_mem_net{f,2}(year_day{y,1},domestic(i))>0)) > length(year_day{y,1})/2
                a1(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},1),inv_mem_net{f,1}(year_day{y,1},ins(i)));
                a2(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},2),inv_mem_net{f,1}(year_day{y,1},ins(i)));
                a3(firm_i,y) = corr(inv_mem_net{f,1}(year_day{y,1},3),inv_mem_net{f,1}(year_day{y,1},ins(i)));
            end
        end
        firm_i = firm_i + 1;
    end
    cc(i,1) = mean(delniz(a1(:),[1,3],0));
    cc(i,2) = mean(delniz(a2(:),[1,3],0));
    cc(i,3) = mean(delniz(a3(:),[1,3],0));
end

%% Inventory correlation with ifd using 'ca', 'cb', 'cc'
% 'ca' is the correlation domestic member(41) with ifd(3)
% 'cb' is the correlation foreign member(21) with ifd(3)
% 'cc' is the correlation institutions(8) with ifd(3)

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
C2 = linspecer(3);

figure1 = figure;
set(gcf,'color','w')
subplot(1,2,1)
scatter(ca(:,1),ca(:,3),mean_vol_member(domestic)/10^9*5,'filled','markerfacecolor',C2(:,1),'MarkerEdgeColor',[0 0 0])
hold on
plot(ca(:,1),ca(:,3),'.','color',[0 0 0])
hold on
plot([-.15 .6],[0 0],'k--')
hold on
plot([0 0],[-.25 .2],'k--')
for i = 1:41
    text(ca(i,1)+.005,ca(i,3),sprintf('%.0f',i),'fontsize', 10)
end
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',14)

xlabel('Corr. with individual','fontsize',16)
ylabel('Corr. with institution','fontsize',16)
title('Domestic member','fontsize',16)
xlim([-.15 .6])
ylim([-.32 .21])

% right upper
a_ru = [11,15,34,32,6,13,30,33,29];
axes('Position',[.65*.5 .65 .3*.5 .25])
scatter(ca(a_ru,1),ca(a_ru,3),mean_vol_member(domestic(a_ru))/10^9*5,'filled','markerfacecolor',C2(:,1),'MarkerEdgeColor',[0 0 0])
hold on
plot(ca(a_ru,1),ca(a_ru,3),'k.')
for i = 1:41
    if sum(a_ru == i) > 0
        text(ca(i,1)+.005,ca(i,3),sprintf('%.0f',i),'fontsize', 10)
    end
end
box on
ylim([.105, .142])
xlim([-.08, .024])

% left lower
a_ll = [35,1,7,27,28,20,24,17,23,14,41,17,39,40];
axes('Position',[.33*.5 .155 .35*.5 .27])
scatter(ca(a_ll,1),ca(a_ll,3),mean_vol_member(domestic(a_ll))/10^9*5,'filled','markerfacecolor',C2(:,1),'MarkerEdgeColor',[0 0 0])
hold on
plot(ca(a_ll,1),ca(a_ll,3),'k.')
for i = 1:41
    if sum(a_ll == i) > 0
        text(ca(i,1)+.005,ca(i,3),sprintf('%.0f',i),'fontsize', 10)
    end
end
box on
xlim([.018, .141])

subplot(1,2,2)
scatter(cb(:,2),cb(:,3),mean_vol_member(foreign)/10^9*5,'filled','markerfacecolor',C2(:,2),'MarkerEdgeColor',[0 0 0])
hold on
plot(cb(:,2),cb(:,3),'.k')
hold on
plot([.12, .295], [0, 0],'k--')
xlim([.12, .295])
ylim([-.08, .01])
for i = 1:21
    if i+41 == 61
        text(cb(i,2)-.01,cb(i,3)+.0012,sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 62
        text(cb(i,2)-.01,cb(i,3)-.0012,sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 57
        text(cb(i,2),cb(i,3)+.002,sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 58
        text(cb(i,2)+.005,cb(i,3),sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 59
        text(cb(i,2)+.005,cb(i,3),sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 60
        text(cb(i,2),cb(i,3)-.003,sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 43
        text(cb(i,2)-.001,cb(i,3)+.003,sprintf('%.0f',i+41),'fontsize',10)
    elseif i+41 == 42
        text(cb(i,2),cb(i,3)+.003,sprintf('%.0f',i+41),'fontsize',10)
    else
        text(cb(i,2)+.005,cb(i,3),sprintf('%.0f',i+41),'fontsize',10)
    end
end
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',14)
xlabel('Corr. with foreigner','fontsize',16)
ylabel('Corr. with institution','fontsize',16)
title('Foreign member','fontsize',16)


axes('Position',[.8 .4 .15/1.3 .3/1.3])
scatter(ca(:,2),ca(:,3),mean_vol_member(domestic)/10^9*5/5,'filled','markerfacecolor',C2(:,1),'MarkerEdgeColor',[0 0 0])
hold on
scatter(cb(:,2),cb(:,3),mean_vol_member(foreign)/10^9*5/5,'filled','markerfacecolor',C2(:,2),'MarkerEdgeColor',[0 0 0])
hold on
plot([0, 0],[-.4, .23],'k--')
hold on
plot([-.45, .4],[0, 0],'k--')
xlim([-.4, .38])
ylim([-.35, .23])
box on




% sort by correlation with inv_mem_net of fore
[~,f1] = sort(cb(:,1),'descend');
[~,d1] = sort(ca(:,1),'descend');
[~,f2] = sort(cb(:,2),'descend');
[~,d2] = sort(ca(:,2),'descend');
[~,f3] = sort(cb(:,3),'descend');
[~,d3] = sort(ca(:,3),'descend');
domestic_sort(:,1) = domestic(d1);
domestic_sort(:,2) = domestic(d2);
domestic_sort(:,3) = domestic(d3);
foreign_sort(:,1) = foreign(f1);
foreign_sort(:,2) = foreign(f2);
foreign_sort(:,3) = foreign(f3);




%% RMT for correlation matrix

num_firm = 1210;

corr_fr = zeros(num_firm,length(year_day));
corr_fr_1 = zeros(num_firm,length(year_day));
corr_fr_2 = zeros(num_firm,length(year_day));
corr_fr_3 = zeros(num_firm,length(year_day));
corr_member_1 = zeros(num_firm,length(year_day));
corr_member_2 = zeros(num_firm,length(year_day));
corr_member_d = zeros(num_firm,length(year_day));
corr_member_f = zeros(num_firm,length(year_day));
herd_adf = cell(num_firm,3);
leval = zeros(num_firm,length(year_day));
vote = zeros(75,1);
vote_count = zeros(75,1);
RMT_valid = zeros(num_firm,length(year_day));
pnh = cell(num_firm,length(year_day)-1);
pnhd = cell(num_firm,length(year_day)-1);
pnh_d = cell(num_firm,length(year_day)-1);
pnhd_d = cell(num_firm,length(year_day)-1);
pnh_f = cell(num_firm,length(year_day)-1);
pnhd_f = cell(num_firm,length(year_day)-1);

ki_not = 0; % Kiwoom is not present
for f = 1:num_firm
    f
    for y = 1:length(year_day)-1
        onoff_trading = find(vol{f,1}(year_day{y,1}) > 0);
        
        active_day = sum(sign(inv_mem_net{f,2}(year_day{y,1},:)),1);
        active_inv = find(active_day > length(onoff_trading)/1.5);
        clear active day onoff_trading
        
        domestic_act = intersect(domestic, active_inv);
        foreign_act = intersect(foreign, active_inv);
        member_act = [domestic_act, foreign_act];
        clear active_inv
        
        if length(domestic_act) + length(foreign_act) > 0
            corr_df = corr(inv_mem_net{f,1}(year_day{y,1},member_act));

            [eig_vec, eig_val] = eig(corr_df);
            [eig_val_sort, ind_val_sort] = sort(on_diagonal(eig_val),'descend');
            largest_eig_val = eig_val_sort(1);
            largest_eig_vec = eig_vec(:,ind_val_sort(1));
            
            Q = size(year_day{y,1},1) / size(member_act,2);
            lambda_min = 1 + 1/Q - 2*sqrt(1/Q);
            lambda_max = 1 + 1/Q + 2*sqrt(1/Q);
            
%             figure
%             subplot(1,2,1)
%             x=linspace(lambda_min,lambda_max,10000); 
%             plot(x,Q/(2*pi)*sqrt((lambda_max-x).*(x-lambda_min))./x)
%             subplot(1,2,2)
%             histfit(eig_val_sort,30)
            
            if sum(member_act==46) > 0
                if largest_eig_vec(find(member_act==46)) <= 0
                    factor = inv_mem_net{f,1}(year_day{y,1},member_act) * largest_eig_vec;
                    vote(member_act) = vote(member_act) + largest_eig_vec;
                elseif largest_eig_vec(find(member_act==46)) > 0
                    factor = inv_mem_net{f,1}(year_day{y,1},member_act) * (-largest_eig_vec);
                    vote(member_act) = vote(member_act) + (-largest_eig_vec);
                end

                vote_count(member_act) = vote_count(member_act) + 1;

                q = length(year_day{y,1}) / length(member_act);
                if largest_eig_val > 1+ 1/q + 2*sqrt(1/q)
                    RMT_valid(f,y) = 1;
                end

                leval(f,y) = largest_eig_val;

                % 'corr_fr' is the correlation coefficient between factor and price_return
                corr_fr(f,y) = corr(factor,price_return{f,1}(year_day{y,1}));

                % 'corr_fr_1' is the correlation coefficient between individual and price_return
                % 'corr_fr_2' is the correlation coefficient between foreigner and price_return
                % 'corr_fr_3' is the correlation coefficient between institution and price_return
                corr_fr_1(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},1),price_return{f,1}(year_day{y,1}));
                corr_fr_2(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},2),price_return{f,1}(year_day{y,1}));
                corr_fr_3(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},3),price_return{f,1}(year_day{y,1}));
                 
                % 'x1' is the sum of direction in daily
                x1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},member_act)),2);
                % 'x1_1' is the sum of positive direction in daily
                x1_1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},member_act))>0,2);
                % 'x2' is the number of member_act in daily
                x2 = sum(sign(inv_mem_net{f,2}(year_day{y,1},member_act)),2);
                % 'x3' is the price return in daily
                x3 = price_return{f,1}(year_day{y,1});
                
                % 'x1_d' is the sum of direction of [domestic member] in daily
                x1_d = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act)),2);
                % 'x1_d_1' is the sum of positive direction of [domestic member] in daily
                x1_d_1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act))>0,2);
                % 'x2_d' is the number of member_act of [domestic member] in daily
                x2_d = sum(sign(inv_mem_net{f,2}(year_day{y,1},domestic_act)),2);
                
                % 'x1_f' is the sum of direction of [foreign member] in daily
                x1_f = sum(sign(inv_mem_net{f,1}(year_day{y,1},foreign_act)),2);
                % 'x1_f_1' is the sum of positive direction of [foreign member] in daily
                x1_f_1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},foreign_act))>0,2);
                % 'x2_f' is the number of member_act of [foreign member] in daily
                x2_f = sum(sign(inv_mem_net{f,2}(year_day{y,1},foreign_act)),2);
                
                
                herd_adf{f,1} = [herd_adf{f,1}; x1./x2];
                herd_adf{f,2} = [herd_adf{f,2}; x1_d./x2_d];
                herd_adf{f,3} = [herd_adf{f,3}; x1_f./x2_f];
                
                % herding using binomial null hypothesis
                % 'pnh' is the existence of herding
                % 'pnhd' is the direction of herding
                for bn = 1:size(year_day{y,1},1)
                    if nchoosek(x2(bn,1),x1_1(bn,1))*(1/2)^(x2(bn,1)) <= 0.05
                        pnh{f,y}(bn,1) = 1;
                        if x1_1(bn,1) < x2(bn,1)/2
                            pnhd{f,y}(bn,1) = -1;
                        else
                            pnhd{f,y}(bn,1) = 1;
                        end
                    else
                        pnh{f,y}(bn,1) = 0;
                        pnhd{f,y}(bn,1) = 0;
                    end
                    
                    if nchoosek(x2_f(bn,1),x1_f_1(bn,1))*(1/2)^(x2_f(bn,1)) <= 0.05
                        pnh_f{f,y}(bn,1) = 1;
                        if x1_f_1(bn,1) < x2_f(bn,1)/2
                            pnhd_f{f,y}(bn,1) = -1;
                        else
                            pnhd_f{f,y}(bn,1) = 1;
                        end
                    else
                        pnh_f{f,y}(bn,1) = 0;
                        pnhd_f{f,y}(bn,1) = 0;
                    end
                    
                    if nchoosek(x2_d(bn,1),x1_d_1(bn,1))*(1/2)^(x2_d(bn,1)) <= 0.05
                        pnh_d{f,y}(bn,1) = 1;
                        if x1_d_1(bn,1) < x2_d(bn,1)/2
                            pnhd_d{f,y}(bn,1) = -1;
                        else
                            pnhd_d{f,y}(bn,1) = 1;
                        end
                    else
                        pnh_d{f,y}(bn,1) = 0;
                        pnhd_d{f,y}(bn,1) = 0;
                    end
                end
            else
                
                ki_not = ki_not + 1;
                
                ki_not_list(ki_not,1) = f;
                ki_not_list(ki_not,2) = y;
            end
            
        end
    end
end


% 'pnh' is the existence of herding
% 'pnhd' is the direction of herding
% 'pnhh' is number of day of herding
% 'pnhhp' is number of day of buy herding
% 'pnhhn' is number of day of sell herding
% 'pnhh_d' is number of day of domestic member herding
% 'pnhh_f' is number of day of foreign member herding
for f = 1:num_firm
    for y = 1:size(year_day,1)-1
        pnhh(f,y) = sum(pnh{f,y}) ./ size(year_day{y,1},1);
        pnhdp(f,y) = sum(pnhd{f,y}>0) ./  size(year_day{y,1},1);
        pnhdn(f,y) = sum(pnhd{f,y}<0) ./  size(year_day{y,1},1);
        
        pnhh_d(f,y) = sum(pnh_d{f,y}) ./ size(year_day{y,1},1);
        pnhh_f(f,y) = sum(pnh_f{f,y}) ./ size(year_day{y,1},1);
    end
end



for f = 1:num_firm
    for y = 1:size(year_day,1)-1
        if ~isempty(pnhd{f,y})
            cdr(f,y) = corr(pnhd{f,y},price_return{f,1}(year_day{y,1}));
        end
        if ~isempty(pnhd_d{f,y})
            cdr_d(f,y) = corr(pnhd_d{f,y},price_return{f,1}(year_day{y,1}));
        end
        if ~isempty(pnhd_f{f,y})
            cdr_f(f,y) = corr(pnhd_f{f,y},price_return{f,1}(year_day{y,1}));
        end
    end
end


for f = 1:10
    pnhh_mean(f,1) = mean(mean(pnhh(1+(f-1)*120:f*120,:)));
    pnhh_d_mean(f,1) = mean(mean(pnhh_d(1+(f-1)*120:f*120,:)));
    pnhh_f_mean(f,1) = mean(mean(pnhh_f(1+(f-1)*120:f*120,:)));
    
    temp = [];
    temp_d = [];
    temp_f = [];
    for i = 1+(f-1)*120:f*120
        temp = [temp; delniz(cdr(i,:),[1,2],0)];
        temp_d = [temp_d; delniz(cdr_d(i,:),[1,2],0)];
        temp_f = [temp_f; delniz(cdr_f(i,:),[1,2],0)];
    end
    cdr_mean(f,1) = mean(temp);
    cdr_d_mean(f,1) = mean(temp_d);
    cdr_f_mean(f,1) = mean(temp_f);
end




figure
set(gcf,'color','w')
subplot(1,2,1)
plot(pnhh_mean,'color',C2(:,1),'linewidth',2)
hold on
plot(pnhh_d_mean,'color',C2(:,2),'linewidth',2)
hold on
plot(pnhh_f_mean,'color',C2(:,3),'linewidth',2)
xlim([1 10])
xlabel('Market capitalization decile')
ylabel('<H>')
legend('All members','Domestic members','Foreign members')

subplot(1,2,2)
plot(cdr_mean,'color',C2(:,1),'linewidth',2)
hold on
plot(cdr_d_mean,'color',C2(:,2),'linewidth',2)
hold on
plot(cdr_f_mean,'color',C2(:,3),'linewidth',2)
hold on
plot([1 10],[0 0],'--k')
xlim([1 10])
xlabel('Market capitalization decile')
ylabel('<DH>')
legend('All members','Domestic members','Foreign members')



% check which direction for kiwoom to move compared with price return
clear aa
for f = 1:121
    for y = 1:11
        aa(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},46),price_return{f,1}(year_day{y,1}));
    end
end
sprintf('direction of Kiwoom: negative %.2f %%',length(find(delniz(aa(:),[1,2,3],0)<0))/length(aa(:))*100)


% timeseries (yearly): factor(individual, foreigner, institution) vs. return
fs = 1;
fe = 121;
clear a b
for i = 1:11
    a(i,1) = mean(delniz(corr_fr(fs:fe,i),[1,2,3],0));
    a(i,2) = mean(delniz(corr_fr_1(fs:fe,i),[1,2,3],0));
    a(i,3) = mean(delniz(corr_fr_2(fs:fe,i),[1,2,3],0));
    a(i,4) = mean(delniz(corr_fr_3(fs:fe,i),[1,2,3],0));
end

figure
plot(a(:,1),'k--')
hold on
plot(a(:,2),'r')
hold on
plot(a(:,3),'g')
hold on
plot(a(:,4),'b')
legend('factor','individual','foreigner','institution')



%% correlation vs. largest eigenvalue

figure
set(gcf,'color','w')

f = 10;
leval_temp = leval((f-1)*121+1 : f*121,:);
corr_fr_temp = corr_fr((f-1)*121+1 : f*121,:);

leval_121 = leval_temp(:);
corr_fr_121 = corr_fr_temp(:);

hold on
scatter(corr_fr_121,leval_121,'b.')
xlabel('corr. factor & return')
ylabel('largest eigen value')

f = 1;
leval_temp = leval((f-1)*121+1 : f*121,:);
corr_fr_temp = corr_fr((f-1)*121+1 : f*121,:);

leval_121 = leval_temp(:);
corr_fr_121 = corr_fr_temp(:);

hold on
scatter(corr_fr_121,leval_121,'r')





figure1 = figure;
set(gcf,'color','w')
subplot(1,2,1)
for ff = 1:10
%     f = 11 - ff;
	f = ff;
    for i = (f-1)*121+1 : f*121
        if f == 1 || f == 10
            scatter(mean(corr_fr(i,1:11)),mean(leval(i,1:11)),'filled','MarkerFaceColor',[f/14 f/14 f/14],'MarkerEdgeColor',[f/14 f/14 f/14])
            hold on
        end
    end
end
text(-0.14,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('Corr. factor and return','Fontsize',12)
ylabel('Largest eigenvalue','Fontsize',12)

subplot(1,2,2)
for f = 1:10
    scatter(mean(mean(corr_fr((f-1)*121+1 : f*121,1:11))),mean(mean(leval((f-1)*121+1 : f*121,1:11))),'MarkerFaceColor',[f/14 f/14 f/14],'MarkerEdgeColor',[f/14 f/14 f/14])
    text(mean(mean(corr_fr((f-1)*121+1 : f*121,1:11)))+.01,mean(mean(leval((f-1)*121+1 : f*121,1:11))),sprintf('%.0f',f))
    hold on
    
    cc(f,1) = mean(mean(corr_fr((f-1)*121+1 : f*121,1:11)));
    cc(f,2) = mean(mean(leval((f-1)*121+1 : f*121,1:11)));
end

hold on
[b, ~, ~, ~,stats] = regress(cc(:,2),[ones(size(cc(:,1))), cc(:,1)]);
pt = plot([-.05 .4],[b(1)+b(2)*(-.05) b(1)+b(2)*(.4)],'k--');
legend(pt(1),sprintf('y=%.2f+%.2fx, R^2=%.2f',b(2),b(1),stats(1)))
legend boxoff
text(-0.14,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('Corr. factor and return','Fontsize',12)
ylabel('Largest eigenvalue','Fontsize',12)
xlim([-.04 .4])

dx=0.08;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 450])



%% Calculating trending 
% Corr of ratio & net_buy with price_return

member_corr_net = cell(10,75);
member_corr_ratio = cell(10,75);

for f = 1:10
    for i = (f-1)*121+1 : f*121
        [f,i]
        for y = 1:11
            for ii = 1:75
                if sum(sign(inv_mem_net{i,2}(year_day{y,1},ii))) > size(year_day{y,1},1)/1.5
                    % for net buy
                    member_corr_net{f,ii} = [member_corr_net{f,ii}; corr(price_return{i,1}(year_day{y,1}), inv_mem_net{i,1}(year_day{y,1},ii))];

                    % for ratio
                    temp = inv_mem_net{i,1}(year_day{y,1},ii)./inv_mem_net{i,2}(year_day{y,1},ii);
                    [res, ind] = delniz(temp,[1,2],1);
                    if length(res) > 3
                        member_corr_ratio{f,ii} = [member_corr_ratio{f,ii}; corr(price_return{i,1}(year_day{y,1}(ind)), res)];
                    end
                end
            end
        end
    end
end


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

k = 3;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
C2 = linspecer(3);


figure1 = figure;
set(gcf,'color','w')
plot(trending(:,1),directionality_member{1,1}(:,k),'-','color',C2(:,1),'linewidth',1.5)
hold on
plot(trending(:,2),directionality_member{2,1}(:,k),'--','color',C2(:,2),'linewidth',1.5)
hold on
plot(trending(:,3),directionality_member{3,1}(:,k),':','color',C2(:,3),'linewidth',1.5)
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




figure
set(gcf,'color','w')
for ii = 1:40
    subplot(4,10,ii)
    plot(trending(:,domestic(ii)),directionality_member{domestic(ii),1}(:,k),'--k','linewidth',1.5)
    hold on
    for i = 1:10
        scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
        hold on
    end
    plot([0 0],[0 1],'k--')
    hold on
    plot([-1 1],[.5 .5],'k--')
    title(inv_mem_name(domestic(ii)))
    xlim([-.5 .2])
    ylim([.4 1])
end

% xlim([-.5 .5])
% ylim([.4 1])



figure
set(gcf,'color','w')

for ii = 1:41
    if ca(ii,3) < -.05 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/5
        subplot(1,3,1)
%         plot(trending(:,domestic(ii)),directionality_member{domestic(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 1:1
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
            hold on
        end
        text(trending(1,domestic(ii)),directionality_member{domestic(ii),1}(1,k),inv_mem_name(domestic(ii)))
    elseif ca(ii,3) > .1 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/5
        subplot(1,3,2)
%         plot(trending(:,domestic(ii)),directionality_member{domestic(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 1:1
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
            hold on
        end
        text(trending(1,domestic(ii)),directionality_member{domestic(ii),1}(1,k),inv_mem_name(domestic(ii)))
    end
end
subplot(1,3,1)
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
xlim([-.5 .5])
ylim([.4 1])
subplot(1,3,2)
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
xlim([-.5 .5])
ylim([.4 1])

subplot(1,3,3)
for ii = 1:21
    if cb(ii,2) > .1 && mean_vol_member(foreign(ii)) > mean(mean_vol_member(foreign))/2
%         plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 1:1
            scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
            hold on
        end
        text(trending(1,foreign(ii)),directionality_member{foreign(ii),1}(1,k),inv_mem_name(foreign(ii)))
    end
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
xlim([-.5 .5])
ylim([.4 1])



%% directionality & trending for ifd

for i = 1:10
    inst_a(i,1) = mean(delniz(dome{i,1},[1,2],0));
    fore_a(i,1) = mean(delniz(fore{i,1},[1,2],0));
    indv_a(i,1) = mean(delniz(indv{i,1},[1,2],0));
end


addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
C2 = linspecer(3);

k = 2;
figure1 = figure;
set(gcf,'color','w')
subplot(1,2,1)
plot(indv_a,lr_indv(:,k),'-','color',C2(:,1),'linewidth',1.5)
hold on
plot(fore_a,lr_fore(:,k),'--','color',C2(:,2),'linewidth',1.5)
hold on
plot(inst_a,lr_inst(:,k),':','color',C2(:,3),'linewidth',1.5)
hold on
plot([0 0],[0 1],'k--')
hold on
plot([-.6 .35],[.5 .5],'k--')
for i = 1:10
    hold on
    scatter(indv_a(i,1),lr_indv(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
    hold on
    scatter(fore_a(i,1),lr_fore(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
    hold on
    scatter(inst_a(i,1),lr_inst(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
end
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.6 .35])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
% lgd1 = legend('Individual','Foreigner','Institution');
lgd1 = legend('Kiwoom','Morgan Stanley','KTB');
lgd1.FontSize = 12;
legend boxoff

subplot(1,2,2)
plot(kl(:,1),'-','color',C2(:,1),'linewidth',1.5)
hold on
plot(kl(:,2),'--','color',C2(:,2),'linewidth',1.5)
hold on
plot(kl(:,3),':','color',C2(:,3),'linewidth',1.5)
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([1 10])
lgd2 = legend('Kiwoom - Morgan Stanley','Kiwoom - KTB','Morgan Stanley - KTB');
lgd2.FontSize = 12;
legend boxoff
xlabel('Market capitalization decile','fontsize',12)
ylabel('Jensen-Shannon divergence','fontsize',12)

dx=0.08;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 450])


%% save network
num_firm = 1216;

for f = 1:num_firm
    for y = 1:length(year_day)
        onoff_trading = find(vol{f,1}(year_day{y,1}) > 0);
        
        active_day = sum(sign(inv_mem_net{f,2}(year_day{y,1},:)),1);
        active_firm = find(active_day > length(onoff_trading)/1.5);
        clear active day onoff_trading
        
        domestic_act = intersect(domestic, active_firm);
        foreign_act = intersect(foreign, active_firm);
        member_act = [domestic_act, foreign_act];
        clear active_firm
        
        if length(domestic_act) + length(foreign_act) > 0
            corr_df = corr(inv_mem_net{f,1}(year_day{y,1},member_act));
            
            % save node
            fid = fopen(sprintf('./save_network/inv_node_%d_%d.csv',f,y),'w');
            fprintf(fid,'id,label,color\n');
            for i = 1:size(domestic_act,2)
                fprintf(fid,'%d,%s,%d\n',domestic_act(i),inv_mem_name{domestic_act(i),1},1);
            end
            for i = 1:size(foreign_act,2)
                fprintf(fid,'%d,%s,%d\n',foreign_act(i),inv_mem_name{foreign_act(i),1},2);
            end
            fclose(fid);
            % save edge
            fid = fopen(sprintf('./save_network/inv_edge_%d_%d.csv',f,y),'w');
            fprintf(fid,'source,target,weight\n');
            for i = 1:size(corr_df,1)-1
                for ii = i+1:size(corr_df,2)
                    fprintf(fid,'%d,%d,%d\n',member_act(i),member_act(ii),corr_df(i,ii));
                end
            end
            fclose(fid);
        end
    end
end


for f = 1:num_firm
    onoff_trading = find(vol{f,1} > 0);

    active_day = sum(sign(inv_mem_net{f,2}),1);
    active_firm = find(active_day > length(onoff_trading)/1.5);
    clear active day onoff_trading

    domestic_act = intersect(domestic, active_firm);
    foreign_act = intersect(foreign, active_firm);
    member_act = [domestic_act, foreign_act];
    clear active_firm

    if length(domestic_act) + length(foreign_act) > 0
        corr_df = corr(inv_mem_net{f,1}(:,member_act));

        % save node
        fid = fopen(sprintf('./save_network/inv_node_%d.csv',f),'w');
        fprintf(fid,'id,label,color\n');
        for i = 1:size(domestic_act,2)
            fprintf(fid,'%d,%s,%d\n',domestic_act(i),inv_mem_name{domestic_act(i),1},1);
        end
        for i = 1:size(foreign_act,2)
            fprintf(fid,'%d,%s,%d\n',foreign_act(i),inv_mem_name{foreign_act(i),1},2);
        end
        fclose(fid);
        % save edge
        fid = fopen(sprintf('./save_network/inv_edge_%d.csv',f),'w');
        fprintf(fid,'source,target,weight\n');
        for i = 1:size(corr_df,1)-1
            for ii = i+1:size(corr_df,2)
                fprintf(fid,'%d,%d,%d\n',member_act(i),member_act(ii),corr_df(i,ii));
            end
        end
        fclose(fid);
    end
end


