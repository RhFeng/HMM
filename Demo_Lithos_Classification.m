clear
close all
% A demo for the lithology classification by Hidden Markov Models
% Input is the inverion result

Lithos_num = 12;
% 1 for CS_non; 2 for MS_non; 3 for MS; 4 for FS_non; 5 for FS
% 6 for VFS_non; 7 for VFS; 8 for SS_non; 9 for SS; 10 for Clay_non
% 11 for Clay; 12 for Coal

% load the inversion results and truth in terms of kappa and M
% the unfiltered truth
% the sampling interval in the vertical direction is 0.4 m

kappa_hc_2D = load('kappa_hc_2D.dat');
M_hc_2D     = load('M_hc_2D.dat');
dz          = 0.4;

% the filtered truth
% the sampling interval in the vertical direction is 5 m
kappa_out   = load('kappa_out.dat');
M_out       = load('M_out.dat');
dz_out      = 5;

% the inversion results
% the sampling interval in the vertical direction is 5 m
kappa_inv   = load('kappa_inv.dat');
M_inv       = load('M_inv.dat');

% load the lithology
% the sampling interval in the vertical direction is 0.4 m
Lithos_unrsa = load('Lithos_unrsa.dat'); %unsampled

% the sampling interval in the vertical direction is 5 m
Lithos_rsa   = load('Lithos_rsa.dat');   %resampled

% choose the logging location for the training of emission functions
CMP_log = [185 1130 1880]; %[185 1130 1880]


% display the inversion results
c_lim_kappa = [min(min(kappa_out)) max(max(kappa_out))];
c_lim_M = [min(min(M_out)) max(max(M_out))];
[Nz_out, N_pan] = size(kappa_out);

figure;
set(gcf,'unit','inches','position',[0.5 0.5 14 6])
subplot('position',[0.07 0.45 0.42 0.3]);
imagesc([1,N_pan],[1,dz_out*Nz_out],kappa_out,c_lim_kappa);
title('True \it\kappa','fontsize',16);
ylabel(' \leftarrow  \itz \rm(m)','fontsize',16);
xlabel('CMP','fontsize',16);
set(gca,'YTick',[0:100:500]);
set(gca,'XTick',[0:400:2000]);
set(gca,'fontsize',20,'linewidth',2);

subplot('position',[0.55 0.45 0.42 0.3]);
imagesc([1,N_pan],[1,dz_out*Nz_out],kappa_inv,c_lim_kappa);
title('Inverted \it\kappa','fontsize',16);
xlabel('CMP','fontsize',16);
set(gca,'YTick',[0:100:500]);
set(gca,'XTick',[0:400:2000]);
set(gca,'fontsize',20,'linewidth',2);

colormap(jet);
h=colorbar('southoutside','position',[0.07 0.20 0.90 0.05]);
t=title(h,'m^2 / N','fontsize',16);
set(get(h,'Title'),'string','m^2 / N','fontsize',16);

figure;
set(gcf,'unit','inches','position',[0.5 0.5 14 6])
subplot('position',[0.07 0.45 0.42 0.3]);
imagesc([1,N_pan],[1,dz_out*Nz_out],M_out,c_lim_M);
title('True \itM','fontsize',16);
ylabel(' \leftarrow  \itz \rm(m)','fontsize',16);
xlabel('CMP','fontsize',16);
set(gca,'YTick',[0:100:500]);
set(gca,'XTick',[0:400:2000]);
set(gca,'fontsize',20,'linewidth',2);

subplot('position',[0.55 0.45 0.42 0.3]);
imagesc([1,N_pan],[1,dz_out*Nz_out],M_inv,c_lim_M);
title('Inverted \itM','fontsize',16);
xlabel('CMP','fontsize',16);
set(gca,'YTick',[0:100:500]);
set(gca,'XTick',[0:400:2000]);
set(gca,'fontsize',20,'linewidth',2);

colormap(jet);
h=colorbar('southoutside','position',[0.07 0.20 0.90 0.05]);
t=title(h,'m^2 / N','fontsize',16);
set(get(h,'Title'),'string','m^2 / N','fontsize',16);

% display the lithology

load mycolor.mat;

figure;
set(gcf,'unit','inches','position',[0.5 0.5 14 5])
subplot('position',[0.07 0.15 0.8 0.75]);
imagesc([1,N_pan],[1,dz_out*Nz_out],Lithos_rsa,[1 12]);
title('Lithology','fontsize',16);
ylabel(' \leftarrow  \itz \rm(m)','fontsize',16);
xlabel('CMP','fontsize',16);
set(gca,'fontsize',20,'linewidth',2);
colormap(mycolor);
labels={'CS\_non','MS\_non','MS','FS\_non','FS','VFS\_non','VFS','SS\_non','SS','Clay\_non','Clay','Coal'};
lcolorbar(labels,'fontweight','bold','fontsize',16);
% plotting the logging location
for i = 1 : length(CMP_log)
    hold on;
    plot([CMP_log(i) CMP_log(i)], [1 dz_out*Nz_out(end)], 'k-.', 'linewidth', 2);
end
hold off;

%% Emission funcion
% build the emission function (Pr(y|x)) based on the selected well logs (CMP_log)
% y: kappa and M
% x: lithology

% firstly select the properties according to CMP_log
% there is tapering effect in the data
% thus we are only choosing the data from 126 to 1175

Prop = [];
factor = 1e11;
for i = 1 : length(CMP_log)
    Prop = [Lithos_unrsa(126:1175,CMP_log(i)) kappa_hc_2D(126:1175,CMP_log(i)).*factor M_hc_2D(126:1175,CMP_log(i)).*factor; Prop];
end

% calculate the mean and covariance for each lithology
% notice that it is arguable that the data have not been normalized since 
% kappa and M have the same unit (factor used before in Prop)

% Gaussian assumption is assumed for the distribution of properties given
% lithology

Lithos_Mean_Cov = zeros(Lithos_num,6);

for i = 1 : Lithos_num
    index = find(Prop(:,1) == i);
    if isempty(index) ~= 1
       temp = (Prop(index,2:3));
       Lithos_Mean_Cov(i,1:2) = mean(temp);
       temp_cov = cov(temp);
       Lithos_Mean_Cov(i,3:6) = [temp_cov(1,:) temp_cov(2,:)];
    end
end


%% The transtion matrix
% the transitional probability matrix will be built based on a selected
% well log location
% and notice that different with emission function
% here the resampled lithology profile will be used (5 m)

CMP_train = 1130; %1130
data_loc = [11 94]; % the tapering effect in the model

% now the lithology index has to be changed in order to only consider the
% lithology at the CMP_train location

Ori_Lithos_Index = unique(Lithos_rsa(data_loc(1):data_loc(2),CMP_train));

Lithos_num_new = length(Ori_Lithos_Index);

[N_log,~] = size(Lithos_rsa(data_loc(1):data_loc(2),CMP_train));

Lithos_train_new = zeros(N_log,1);

for i = 1 : N_log
    Lithos_train_new(i,1) = find(Lithos_rsa(data_loc(1)+i-1,CMP_train) == Ori_Lithos_Index);
end

% vertical transition matrix

PV = zeros(Lithos_num_new,Lithos_num_new);

[~, PV] = tp(Lithos_num_new,Lithos_train_new,1);

% find zero entries in the matrix 
% assign eta to exclude abrupt transitions

eta = 0.0001;

for i = 1 : Lithos_num_new
    m = 0;
    for j = 1 : Lithos_num_new
        if PV(i,j) == 0
            PV(i,j) = eta;
            m = m + 1;
        end
    end
    
    a = find(PV(i,:) == max(PV(i,:)));
    
    PV(i,a(1)) = PV(i,a(1)) - m * eta;
end

%% Training based on the inversion results
input_data = [kappa_inv(data_loc(1):data_loc(2),CMP_train) M_inv(data_loc(1):data_loc(2),CMP_train)]';
input_data = input_data * factor;

O = 2; % number of input features (kappa and M)
Q = length(Ori_Lithos_Index); % number of hidden variables (lithologies)
M = length(Ori_Lithos_Index); % number of mixture components
mixmat0 = ones(length(Ori_Lithos_Index),1);

max_iter=100; % number of the training iterations
adj_Prior=0;  % adjust prior (1) or not (0)
adj_trans=0;  % adjust transition matrix (1) or not (0)
adj_mu=1;     % adjust means (1) or not (0)
adj_sigma=1;  % adjust covariance (1) or not (0)
adj_mix=0;    % adjust mixture component (1) or not (0)


for i=1:length(Ori_Lithos_Index)
    temp_mu0(:,i)=Lithos_Mean_Cov(Ori_Lithos_Index(i),1:2)';
    temp_sigma0(:,:,i)=[Lithos_Mean_Cov(Ori_Lithos_Index(i),3:4);Lithos_Mean_Cov(Ori_Lithos_Index(i),5:6)];
end

P_int=p2stat(PV); % initial lithology distribution



[LL, prior1, PV1, temp_mu1, temp_sigma1,mixmat1] = mhmm_em(input_data, P_int, PV, temp_mu0, temp_sigma0,mixmat0,...
                                                         'max_iter', max_iter,'adj_Prior',adj_Prior,'adj_trans',adj_trans,...
                                                         'adj_mu',adj_mu,'adj_sigma',adj_sigma,'adj_mix',adj_mix,'cov_type','full');

for i = 1 : length(Ori_Lithos_Index)
    for j = 1 : O
        for k = 1 : O
            temp_sigma1(k,j,i) = roundn(temp_sigma1(k,j,i),-5);
        end
    end
end                                                    
                                                     
mu0=zeros(O,Lithos_num);
mu1=zeros(O,Lithos_num);
sigma0=zeros(O,O,Lithos_num);
sigma1=zeros(O,O,Lithos_num);                                                     
                                                     
for i = 1 : length(Ori_Lithos_Index)
    k = Ori_Lithos_Index(i);
    mu0(:,k) = temp_mu0(:,i);
    mu1(:,k) = temp_mu1(:,i);
    
    sigma0(:,:,k) = temp_sigma0(:,:,i);
    sigma1(:,:,k) = temp_sigma1(:,:,i);
end 

legend_string = ['CS\_non  ';'MS\_non  ';'MS       ';'FS\_non  ';'FS       ';'VFS\_non ';'VFS      ';...
    'SS\_non  ';'SS       ';'Clay\_non';'Clay     ';'Coal     '];

cell_string = cellstr(legend_string);
chr = char(cell_string);

legstr = [];

for i = 1 : length(Ori_Lithos_Index)
    ii = Ori_Lithos_Index(i);
    legstr{i} = chr(ii,:);
end

clear h;
cmap = mycolor;

% plotting the probability distribution of 12 lithologies
figure;
set(gcf,'unit','inches','position',[0.5 0.5 6 6])
subplot('position',[0.1 0.12 0.8 0.8]);
for i = 1 : Lithos_num
    
    Z=calcov2(Lithos_Mean_Cov(i,1:2)',[Lithos_Mean_Cov(i,3:4);Lithos_Mean_Cov(i,5:6)]);
    h(i)=plot(Z(1,:),Z(2,:),'color',cmap(i,:),'LineWidth',3);

    hold on;
    title('Emission Function of 12 Lithologies','fontsize',22);
    ylabel('\itM \rm(\times10^{-11})','fontsize',20);
    xlabel('\it\kappa \rm(\times10^{-11})','fontsize',20);
    grid on;
    xlim([2 8]);
    ylim([0,4]);
    set(gca,'fontsize',16);

end
legend(h,chr);
hold off;
set(gca,'fontsize',16,'GridLineStyle','-.','linewidth',2);

figure;
set(gcf,'unit','inches','position',[0.5 0.5 12 6])
subplot('position',[0.08 0.12 0.4 0.8]);
for i=1:length(Ori_Lithos_Index)
    logg=Ori_Lithos_Index(i);
    
    Z=calcov2(mu0(:,logg),reshape(sigma0(:,:,logg),O,O));

    h(logg)=plot(Z(1,:),Z(2,:),'color',cmap(logg,:),'LineWidth',3);

    hold on;
    title('Before training','fontsize',22);
    ylabel('\itM \rm(\times10^{-11})','fontsize',20);
    xlabel('\it\kappa \rm(\times10^{-11})','fontsize',20);
    grid on;
    xlim([0 12]);
    ylim([0,8]);
    set(gca,'fontsize',16);

end
legend(h(Ori_Lithos_Index),chr(Ori_Lithos_Index,:));
hold off;
set(gca,'fontsize',16,'GridLineStyle','-.','linewidth',2);

subplot('position',[0.55 0.12 0.4 0.8]);
for i=1:length(Ori_Lithos_Index)
    logg=Ori_Lithos_Index(i);
    
    Z=calcov2(mu1(:,logg),reshape(sigma1(:,:,logg),O,O));

    h(logg)=plot(Z(1,:),Z(2,:),'color',cmap(logg,:),'LineWidth',3);

    hold on;
%     legend(h(logg),chr(logg,:));
    title('After training','fontsize',22);
    ylabel('\itM \rm(\times10^{-11})','fontsize',20);
    xlabel('\it\kappa \rm(\times10^{-11})','fontsize',20);
    grid on;
    xlim([0 12]);
    ylim([0,8]);
    set(gca,'fontsize',16);

end
legend(h(Ori_Lithos_Index),chr(Ori_Lithos_Index,:));
hold off;
set(gca,'fontsize',16,'GridLineStyle','-.','linewidth',2);

%% Now we can do the prediction based on the trained HMM
% first the prediction is performed based on training location (CMP_train)
% then the trained HMM can be used to any close location to CMP_train such
% as CMP = 1140

input_index = 1;

while input_index == 1
    input_CMP = input('Please select a CMP location (it should be between 1130 and 1140): ');
    if input_CMP >= 1130 && input_CMP <= 1140
        input_index = 0;
    else
        disp('Error: the CMP number should be between 1130 and 1140, please re-select!');
    end
end

CMP_predict = input_CMP;


disp(['CMP of Lithology Prediction = ' int2str(CMP_predict)]);

disp('Using the inversion!');

kappa_input = kappa_inv(data_loc(1):data_loc(2),CMP_predict).*factor;
M_input     = M_inv(data_loc(1):data_loc(2),CMP_predict).*factor;

depth = linspace(data_loc(1),data_loc(2),data_loc(2)-data_loc(1)+1) .* dz_out;

figure;
set(gcf,'unit','inches','position',[0.5 0.5 6 8])
subplot('position',[0.15 0.12 0.3 0.8]);
plot(kappa_out(data_loc(1):data_loc(2),CMP_predict),depth,'r','linewidth',3); % truth
ylim([data_loc(1)*dz_out data_loc(2)*dz_out]);
set(gca,'YDir','reverse');
hold on;
plot(kappa_inv(data_loc(1):data_loc(2),CMP_predict),depth,'b','linewidth',3); % inversion
ylim([data_loc(1)*dz_out data_loc(2)*dz_out]);
ylabel('\leftarrow  \itz \rm(m)','fontsize',16);
set(gca,'YDir','reverse');
set(gca,'fontsize',16);
% xlabel('\it\kappa \rm(m^2/N)','fontsize',16);
xlabel('\it\kappa','fontsize',18);
grid on;
set(gca,'fontsize',16,'GridLineStyle','-.','linewidth',2);

set(gcf,'unit','inches','position',[0.5 0.5 6 8])
subplot('position',[0.6 0.12 0.3 0.8]);
plot(M_out(data_loc(1):data_loc(2),CMP_predict),depth,'r','linewidth',3); % truth
ylim([data_loc(1)*dz_out data_loc(2)*dz_out]);
set(gca,'YDir','reverse');
hold on;
plot(M_inv(data_loc(1):data_loc(2),CMP_predict),depth,'b','linewidth',3); % inversion
ylim([data_loc(1)*dz_out data_loc(2)*dz_out]);
% ylabel('\leftarrow  \itz \rm(m)','fontsize',16);
set(gca,'YDir','reverse');
set(gca,'fontsize',16);
% xlabel('\it\kappa \rm(m^2/N)','fontsize',16);
xlabel('\itM','fontsize',18);
grid on;
set(gca,'fontsize',16,'GridLineStyle','-.','linewidth',2);

% now the emission function is being calculated based on the selected
% inversion results
Emiss = zeros(length(M_input),length(Ori_Lithos_Index));

for i = 1 : length(M_input)
    for j = 1 : length(Ori_Lithos_Index)
        Emiss(i,j) = mvnpdf([kappa_input(i) M_input(i)],...
        mu1(:,Ori_Lithos_Index(j))',...
        reshape(sigma1(:,:,Ori_Lithos_Index(j)),O,O));
    end
end

[alpha,beta,gamma,~] = fwdback(prior1,PV1,Emiss');

% prediction with Markov prior
Lithos_predict = zeros(1,length(M_input));

% prediction without Markov prior
Lithos_emiss   = zeros(1,length(M_input));

for i = 1 : length(M_input)
    Lithos_predict(1,i) = find(gamma(:,i) == max(gamma(:,i)));
    Lithos_predict(1,i) = Ori_Lithos_Index(Lithos_predict(1,i));
    
    Lithos_emiss(1,i) = find(Emiss(i,:) == max(Emiss(i,:)));
    Lithos_emiss(1,i) = Ori_Lithos_Index(Lithos_emiss(1,i));
end

figure;
set(gcf,'unit','inches','position',[0.5 0.5 10 8])
subplot('position',[0.1 0.12 0.2 0.8]);
imagesc([0,1],depth,Lithos_rsa(data_loc(1):data_loc(2),CMP_predict),[1 12]);
title('Truth','fontsize',16);
set(gca,'XTickLabel',' ');
ylabel('\leftarrow  \itz \rm(m)','fontsize',16)
set(gca,'fontsize',16,'linewidth',2);

subplot('position',[0.38 0.12 0.2 0.8]);
imagesc([0,1],depth,Lithos_emiss',[1 12]);
title('Prediction (no Markov)','fontsize',16);
set(gca,'XTickLabel',' ');
set(gca,'fontsize',16,'linewidth',2);

subplot('position',[0.66 0.12 0.2 0.8]);
imagesc([0,1],depth,Lithos_predict',[1 12]);
title('Prediction (HMM)','fontsize',16);
set(gca,'XTickLabel',' ');
set(gca,'fontsize',16,'linewidth',2);

colormap(cmap);
labels={'CS\_non','MS\_non','MS','FS\_non','FS','VFS\_non','VFS','SS\_non','SS','Clay\_non','Clay','Coal'};
% labels={' ','MS_non','MS','FS_non','FS','VFS_non','VFS',' ','SS','Clay_non','Clay',' '};
lcolorbar(labels,'fontweight','bold','fontsize',16);
set(gca,'XTickLabel',' ');
set(gca,'fontsize',16);

% confusion matrix
% prediction without Markov
class_emiss = confusion_matrix(12,Lithos_rsa(data_loc(1):data_loc(2),CMP_predict),Lithos_emiss',1);
MCC_emiss      = MCC_coeff(class_emiss)

% confusion matrix
% prediction with Markov
class_predict = confusion_matrix(12,Lithos_rsa(data_loc(1):data_loc(2),CMP_predict),Lithos_predict',1);
MCC_predict   = MCC_coeff(class_predict)

% Conclusion
% it can be seen that with an incorporation of Markov prior, the prediction
% is getting better, as well as reflected by an increase of the correlation
% coefficient (MCC).

