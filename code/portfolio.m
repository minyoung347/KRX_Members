
%% Portfolio Asymmetry

ind_target = ind_kosdaq;
target_return = kosdaq_return;
 
% domestic_ind = find(ca(:,3)<0 & ca(:,1)>0);
% domestic_ins = find(ca(:,3)>=0 & ca(:,1)<0);
domestic_ind = find(ca(:,3)<-0.01);
domestic_ins = find(ca(:,3)>0.01);

asym_domestic = zeros(2801,41);
asym_foreign = zeros(2801,21);
asym_domestic_ind = zeros(2801,length(domestic_ind));
asym_domestic_ins = zeros(2801,length(domestic_ins));

for i = 1:41
    for ii = 1:2801
        asym_domestic(ii,i) = length(find(sign(firm_net{domestic(i),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{domestic(i),1}(ii,ind_target))~=0));
    end
end

for i = 1:21
    for ii = 1:2801
        asym_foreign(ii,i) = length(find(sign(firm_net{foreign(i),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{foreign(i),1}(ii,ind_target))~=0));
    end
end

for i = 1:length(domestic_ind)
    for ii = 1:2801
        asym_domestic_ind(ii,i) = length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,ind_target))~=0));
    end
end

for i = 1:length(domestic_ins)
    for ii = 1:2801
        asym_domestic_ins(ii,i) = length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,ind_target))~=0));
    end
end
        


%%

th=0.00;

for i = 1:50

    asym_d1(i,1) = length(find(delniz(asym_domestic_ind(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ind(:),[1,2],0));
    asym_d2(i,1) = length(find(delniz(asym_domestic_ins(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ins(:),[1,2],0));
    asym_f(i,1) = length(find(delniz(asym_foreign(:),[1,2],0) < 0.5 -th | delniz(asym_foreign(:),[1,2],0) > 0.5 +th))/length(delniz(asym_foreign(:),[1,2],0));
    asym_d(i,1) = length(find(delniz(asym_domestic(:),[1,2],0) < 0.5 -th | delniz(asym_domestic(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic(:),[1,2],0));

    th = th + 0.01;

end

figure
set(gcf,'color','w')
plot(asym_f,'k')
hold on
plot(asym_d2)
hold on
plot(asym_d1)
hold on
plot(asym_d)
% set(gca,'yscale','log')
legend('f','ins','ind','d')
xticks([10 20 30 40 50])
xticklabels({'0.1','0.2','0.3','0.4','0.5'})
legend('Foreign','Domestic(Institutions)','Domestic(Individuals)','Domestic')
xlabel('threshold')
ylabel('asymmetry')

%%

th = 0.00;
cr = zeros(50,3);

for t = 1:50

    for i = 1:2801
        d(i,1) = length(find(delniz(asym_domestic(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic(i,:),[1,2],0));
        di(i,1) = length(find(delniz(asym_domestic_ind(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ind(i,:),[1,2],0));
        ds(i,1) = length(find(delniz(asym_domestic_ins(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ins(i,:),[1,2],0));
        fo(i,1) = length(find(delniz(asym_foreign(i,:),[1,2],0) < 0.5 -th | delniz(asym_foreign(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_foreign(i,:),[1,2],0));
    end

    cr(t,1) = corr(di,abs(target_return));
    cr(t,2) = corr(ds,abs(target_return));
    cr(t,3) = corr(d,abs(target_return));
    cr(t,4) = corr(fo,abs(target_return));

    th = th + 0.01;
    
    aaa = fitlm([di ds fo],abs(target_return));
    bbb(t,1) = aaa.Rsquared.Adjusted;
    
end

% figure
% set(gcf,'color','w')
% plot(cr,'linewidth',1.5)
% xticks([10 20 30 40 50])
% xticklabels({'0.1','0.2','0.3','0.4','0.5'})
% legend('Domestic(Individuals)','Domestic(Institutions)','Domestic','Foreign')
% xlabel('threshold')
% ylabel('Correlation with KOSPI return')

[ta tb] = max(bbb);
th = (tb-1)/100;
for i = 1:2801
    d(i,1) = length(find(delniz(asym_domestic(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic(i,:),[1,2],0));
    di(i,1) = length(find(delniz(asym_domestic_ind(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ind(i,:),[1,2],0));
    ds(i,1) = length(find(delniz(asym_domestic_ins(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ins(i,:),[1,2],0));
    fo(i,1) = length(find(delniz(asym_foreign(i,:),[1,2],0) < 0.5 -th | delniz(asym_foreign(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_foreign(i,:),[1,2],0));
end

fitlm([di ds fo],abs(target_return))






%% Portfolio Asymmetry (10 deciles)

asym_domestic = cell(10,1);
asym_foreign = cell(10,1);
asym_domestic_ind = cell(10,1);
asym_domestic_ins = cell(10,1);

for f = 1:10
    firm_start = (f-1)*121+1;
    firm_end = f*121;
    
    domestic_ind = find(ca(:,3)<0 & ca(:,1)>0);
    domestic_ins = find(ca(:,3)>=0 & ca(:,1)<0);
    
    asym_domestic{f,1} = zeros(2801,41);
    asym_foreign{f,1} = zeros(2801,21);
    asym_domestic_ind{f,1} = zeros(2801,length(domestic_ind));
    asym_domestic_ins{f,1} = zeros(2801,length(domestic_ins));

    for i = 1:41
        for ii = 1:2801
            asym_domestic{f,1}(ii,i) = length(find(sign(firm_net{domestic(i),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{domestic(i),1}(ii,firm_start:firm_end))~=0));
        end
    end

    for i = 1:21
        for ii = 1:2801
            asym_foreign{f,1}(ii,i) = length(find(sign(firm_net{foreign(i),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{foreign(i),1}(ii,firm_start:firm_end))~=0));
        end
    end

    for i = 1:length(domestic_ind)
        for ii = 1:2801
            asym_domestic_ind{f,1}(ii,i) = length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,firm_start:firm_end))~=0));
        end
    end

    for i = 1:length(domestic_ins)
        for ii = 1:2801
            asym_domestic_ins{f,1}(ii,i) = length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,firm_start:firm_end))~=0));
        end
    end
        
end

%%

figure
set(gcf,'color','w')

for f = 1:10

    th=0.00;

    for i = 1:50

        asym_d1(i,1) = length(find(delniz(asym_domestic_ind{f,1}(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind{f,1}(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ind{f,1}(:),[1,2],0));
        asym_d2(i,1) = length(find(delniz(asym_domestic_ins{f,1}(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins{f,1}(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ins{f,1}(:),[1,2],0));
        asym_f(i,1) = length(find(delniz(asym_foreign{f,1}(:),[1,2],0) < 0.5 -th | delniz(asym_foreign{f,1}(:),[1,2],0) > 0.5 +th))/length(delniz(asym_foreign{f,1}(:),[1,2],0));

        th = th + 0.01;

    end

    subplot(2,5,f)
    plot(asym_f)
    hold on
    plot(asym_d2)
    hold on
    plot(asym_d1)
    set(gca,'yscale','log')
    legend('f','ins','ind')
    
end


%%

th = 0.10;
f = 1;

for i = 1:2801
    di(i,1) = length(find(delniz(asym_domestic_ind{f,1}(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind{f,1}(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ind{f,1}(i,:),[1,2],0));
    ds(i,1) = length(find(delniz(asym_domestic_ins{f,1}(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins{f,1}(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ins{f,1}(i,:),[1,2],0));
    fo(i,1) = length(find(delniz(asym_foreign{f,1}(i,:),[1,2],0) < 0.5 -th | delniz(asym_foreign{f,1}(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_foreign{f,1}(i,:),[1,2],0));
end


figure
set(gcf,'color','w')
subplot(1,3,1)
scatter(di,abs(kospi_return),'.')
title(sprintf('%.3f',corr(di,abs(kospi_return))))
set(gca,'yscale','log')
subplot(1,3,2)
scatter(ds,abs(kospi_return),'.')
title(sprintf('%.3f',corr(ds,abs(kospi_return))))
set(gca,'yscale','log')
subplot(1,3,3)
scatter(fo,abs(kospi_return),'.')
title(sprintf('%.3f',corr(fo,abs(kospi_return))))
set(gca,'yscale','log')


figure
set(gcf,'color','w')
scatter(di,abs(kospi_return),'.r')
hold on
scatter(ds,abs(kospi_return),'.g')
hold on
scatter(fo,abs(kospi_return),'.k')
set(gca,'yscale','log')







