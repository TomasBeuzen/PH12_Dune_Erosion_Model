clear
close all force
clc
format bank
format compact

addpath(genpath('./Functions'))

load('june2016.mat')
load('calib_rB05_B0m.mat');

% Select profiles that did erode
idx_selection = find(~calib.idx_noerosion);
% idx_selection = 687;
% idx_selection = find(strcmp([june2016.profiles.site],'BOOM'));
% idx_selection = idx_selection(~ismember(idx_selection,find(calib.idx_noerosion)));
Np = length(idx_selection);

idx_nochange = zeros(Np,1);
Cs_dx = nan(Np,1);
Cs_dV = nan(Np,1);
delta_Cs = nan(Np,1);
idx_best_dx = nan(Np,1);
idx_best_dV = nan(Np,1);
BSS_dx = nan(Np,1);
BSS_dV = nan(Np,1);
% out = f_PH12model_single_v3(june2016,276,0.005,0.5,1)

for i = 1:Np
    
    j = idx_selection(i);
    
    if max(abs(calib.dx(j,:))) - min(abs(calib.dx(j,:))) <= 0.5
       idx_nochange(i) = 1; 
    else   
        temp = find(abs(calib.dx(j,:)) == min(abs(calib.dx(j,:))));
        idx_best_dx(i) = temp(ceil(length(temp)/2));
        Cs_dx(i) = calib.Cs_values(idx_best_dx(i));
        BSS_dx(i) = calib.BSS(j,idx_best_dx(i));
        
        temp = find(abs(calib.dV(j,:)) == min(abs(calib.dV(j,:))));
        idx_best_dV(i) = temp(ceil(length(temp)/2));
        Cs_dV(i) = calib.Cs_values(idx_best_dV(i));
        BSS_dV(i) = calib.dV_norm(j,idx_best_dV(i));
        
        delta_Cs(i) = (log10(Cs_dx(i)) - log10(Cs_dV(i)));
        
%         skip_runs = logical(calib.idx_skip(j,:));        
%         figure(i)
%         h1(i) = subplot(2,1,1);
%         hold on
%         grid on
%         title(sprintf('%s \n %s',[june2016.profiles(j).site{1} ' ' num2str(j)], '\beta_t = 0.5 \beta_0'))
%         yyaxis left
%         hold on
%         plot(calib.Cs_values,calib.BSS(j,:),'b-')
%         plot(calib.Cs_values(skip_runs),calib.BSS(j,skip_runs),'g.')
%         ylim([0 1])
%         plot([Cs_dx(i),Cs_dx(i)],[0 1],'k:')
%         text(3e-5,0.8,['log(\DeltaCs) = ' num2str(abs(delta_Cs(i)),2)],'FontWeight','normal')
%         ylabel('\Deltax')
%         yyaxis right
%         hold on
%         plot(calib.Cs_values,calib.dV_norm(j,:),'r-')
%         plot(calib.Cs_values(skip_runs),calib.dV_norm(j,skip_runs),'g.')
%         ylim([0 1])
%         plot([Cs_dV(i),Cs_dV(i)],[0 1],'k:')
%         ylabel('\DeltaV')
%         set(gca,'XScale','log')
%         set(gca,'XTickLabel',[])
%         pause()
%         close(gcf)

%         out = f_PH12model_single_v3(june2016,j,Cs_dx(i),0.5,1)
%         print(gcf,'out_dx_05','-djpeg','-r300')
%         out = f_PH12model_single_v3(june2016,j,Cs_dV(i),0.5,1)
%         print(gcf,'out_dV_05','-djpeg','-r300')


    end   
end

calib.BSSdx_thresh = 0.8;
idx_keep = ~((abs(delta_Cs) >= 1) | logical(idx_nochange) | BSS_dx < calib.BSSdx_thresh);
calib.idx_keep = idx_keep;
calib.Cs_dx = Cs_dx;
calib.Cs_dV = Cs_dV;
calib.BSS_dV = BSS_dV;
calib.BSS_dx = BSS_dx;

% idx_keep = calib.idx_keep;
sum(calib.idx_keep)
% figure
% h3 = subplot(2,1,1);
% hold on; grid on; box on;
% hHist = histogram(abs(delta_Cs(idx_keep)),'BinWidth',0.05);
% hHist.FaceColor = [.5 .5 .5];
% hHist.LineWidth = 0.5; 
% title('log\DeltaCs for rB = 0.5')
% xlabel('log\DeltaCs');ylabel('counts')
% ylim([0 150])
% print(gcf,'logDeltaCs_rB05','-djpeg','-r300')
%%
% clear calib
% load('calib2_final.mat');
% 
% for i = 1:Np
%     
%     j = idx_selection(i);
%     
%     if max(abs(calib.dx(j,:))) - min(abs(calib.dx(j,:))) <= 0.5
%        idx_nochange(i) = 1; 
%     else   
%         temp = find(abs(calib.dx(j,:)) == min(abs(calib.dx(j,:))));
%         idx_best_dx(i) = temp(ceil(length(temp)/2));
%         Cs_dx(i) = calib.Cs_values(idx_best_dx(i));
%         BSS_dx(i) = calib.BSS(j,idx_best_dx(i));
%         
%         idx_best_dV(i) = find(abs(calib.dV(j,:)) == min(abs(calib.dV(j,:))),1,'first');
%         Cs_dV(i) = calib.Cs_values(idx_best_dV(i));
%         BSS_dV(i) = calib.dV_norm(j,idx_best_dV(i));
%         
%         delta_Cs(i) = (log10(Cs_dx(i)) - log10(Cs_dV(i)));
%         
% %         skip_runs = logical(calib.idx_skip(j,:));      
% %         figure(i)
% %         h2(i) = subplot(2,1,2);
% %         hold on
% %         grid on
% %         title('\beta_t = 1 \beta_0')
% %         yyaxis left
% %         hold on
% %         plot(calib.Cs_values,calib.BSS(j,:),'b-')
% %         plot(calib.Cs_values(skip_runs),calib.BSS(j,skip_runs),'g.')
% %         ylim([0 1])
% %         plot([Cs_dx(i),Cs_dx(i)],[0 1],'k:')
% %         text(3e-5,0.8,['log(\DeltaCs) = ' num2str(abs(delta_Cs(i)),2)],'FontWeight','normal')
% %         ylabel('\Deltax')
% %         yyaxis right
% %         hold on
% %         plot(calib.Cs_values,calib.dV_norm(j,:),'r-')
% %         plot(calib.Cs_values(skip_runs),calib.dV_norm(j,skip_runs),'g.')
% %         ylim([0 1])
% %         plot([Cs_dV(i),Cs_dV(i)],[0 1],'k:')
% %         ylabel('\DeltaV')
% %         set(gca,'XScale','log')
% %         xlim([1e-5 1])
% %         linkaxes([h1(i) h2(i)],'x')
% %         print(figure(i),fullfile('.\figures\Cs_calib_curves\B0m',[june2016.profiles(j).site{1} '_' num2str(j)]),'-djpeg','-r300')
% %         pause()
% %         close(gcf)
%         
% %         out = f_PH12model_single_v3(june2016,j,Cs_dx(i),1,1)
% %         print(gcf,'out_dx_1','-djpeg','-r300')
% %         out = f_PH12model_single_v3(june2016,j,Cs_dV(i),1,1)
% %         print(gcf,'out_dV_1','-djpeg','-r300')
%         
% 
%     end
%         
%         
% end
% 
% % idx_keep = ~((abs(delta_Cs) >= 1) | logical(idx_nochange) | BSS_dx < 0.5)
% % sum(idx_keep)
% % h4 = subplot(2,1,2);
% % hold on; grid on; box on;
% % hHist = histogram(abs(delta_Cs(idx_keep)),'BinWidth',0.05);
% % hHist.FaceColor = [.5 .5 .5];
% % hHist.LineWidth = 0.5; 
% % title('log\DeltaCs for rB = 1')
% % xlabel('log\DeltaCs');ylabel('counts')
% % ylim([0 150])
% % print(gcf,'logDeltaCs_rB1','-djpeg','-r300')
%%
M = calib.M(idx_selection(idx_keep),:);

M(:,13) = Cs_dx(idx_keep);
M(:,14) = Cs_dV(idx_keep);
M(:,15) = log10(Cs_dx(idx_keep));
M(:,16) = log10(Cs_dV(idx_keep));

fields = {'id', 'zb1', 'maxTWL', 'maxOWL', 'profile_orientation',...
        'wave_direction', 'B0', 'B0m', 'dV1', 'dV2', 'maxTWL-zb1', 'maxSWL-zb1',...
        'Cs_dx','Cs_dV','log10Cs_dx','log10Cs_dV'};
% % Save into csv file
% cHeader = fields;
% textHeader = strjoin(cHeader, ',');
% fid = fopen('calib4_final.csv','w'); 
% fprintf(fid,'%s\n',textHeader);
% dlmwrite('calib4_final.csv',M,'-append');
% fclose(fid);

% Make some figures
for i = 11 % 1:12
  
    
    Cs = M(:,13);
    idx = Cs < 4e-1;
    Cs = Cs(idx);
    variable = M(idx,i);
    log10Cs = M(idx,15);
    
    figure
    hold on; grid on; box on;
    set(gca,'XScale','log')
    
    [r] = corrcoef(log10Cs,variable);
    [p] = polyfit(log10Cs,variable,1);
    yfit = polyval(p,log10Cs);
    yresid = variable - yfit; SSresid = sum(yresid.^2); SStotal = (length(variable)-1) * var(variable);
    rsq = 1 - SSresid/SStotal; 
    str = sprintf(' y = %.4f log(x) + %.4f \n R^2 = %.2f   ,   n = %d',p(1),p(2),rsq,length(variable));
    plot(Cs,variable,'k.','displayname','data')
    plot(Cs,yfit,'r-','displayname',str)
    
%     title([fields{i} ' , \rho = ' num2str(r(1,2),2)])
    title(sprintf('Correlation plot, rho = %.2f \n best dx  ,  rB = 1  ,  R2(B0)', r(1,2)))
    xlabel('Cs values');ylabel(fields{i})
    legend('show')
    
%     print(gcf,'corr_dx_rb1_B0','-djpeg','-r300')
end
calib.regression.dx.p1 = p(1);
calib.regression.dx.p2 = p(2);
% Make some figures
for i = 11 % 1:12
  
    
    Cs = M(:,14);
    idx = Cs < 4e-1;
    Cs = Cs(idx);
    variable = M(idx,i);
    log10Cs = M(idx,16);
    
    figure
    hold on; grid on; box on;
    set(gca,'XScale','log')
    
    [r] = corrcoef(log10Cs,variable);
    [p] = polyfit(log10Cs,variable,1);
    yfit = polyval(p,log10Cs);
    yresid = variable - yfit; SSresid = sum(yresid.^2); SStotal = (length(variable)-1) * var(variable);
    rsq = 1 - SSresid/SStotal; 
    str = sprintf(' y = %.4f log(x) + %.4f \n R^2 = %.2f   ,   n = %d',p(1),p(2),rsq,length(variable));
    plot(Cs,variable,'k.','displayname','data')
    plot(Cs,yfit,'r-','displayname',str)
    
%     title([fields{i} ' , \rho = ' num2str(r(1,2),2)])
    title(sprintf('Correlation plot, rho = %.2f \n best dV  ,  rB = 0.5  ,  R2(B0m)', r(1,2)))
    xlabel('Cs values');ylabel(fields{i})
    legend('show')
    
%     print(gcf,'corr_dV_rb1_B0','-djpeg','-r300')
end
calib.regression.dV.p1 = p(1);
calib.regression.dV.p2 = p(2);

figure 
histogram(log10Cs)

