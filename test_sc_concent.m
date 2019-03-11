addpath('common')

%B = [1,2,3,4,5,60;
 %    2,3,4,5,6,1;
  %   3,4,5,6,1,2;
   %  4,5,6,1,2,3];
 
%  B = [10,5,1;
%       11,3,1;
%       12,1,2];

B = 0.5*ones(3) + 2.5*eye(3);
B = [B,0.5*ones(3,1)];


 
[K,L] = size(B);

%n0_vec = [50 100 200 500 1000];
%n0_len = length(n0_vec);
n0_len =  12;
n0_vec = 20*round(logspace(log10(30), log10(2000), n0_len));


%%
T = 3;
%RegVec = [0.5,1,1.2,1.4,inf];
RegVec = [1,1.2,1.4,inf];
nReg = length(RegVec);

acc_all = zeros(n0_len, nReg, 3, T);
nmi_all = zeros(n0_len, nReg, 3, T);
%nmi_all = zeros(n0_len, nReg, T);
dt_all = zeros(n0_len, nReg, T);

prefactor = 1;
al = .5;
beta = 1;
%%
%RegVec = [.5,0.75,1,inf];
for t = 1:T  % This is parallel for; use with caution
    dt = zeros(n0_len,nReg);
    %nmi = zeros(n0_len,nReg,2);
    nmi = zeros(n0_len,nReg,3);
    acc = zeros(n0_len,nReg,3);
    
    for j = 1:n0_len
        fprintf('--- t = %3d, j = %3d ---\n',t,j)
        n0 = n0_vec(j);
        m = L*n0;
        n = K*n0;
        compute_overll_acc = @(y,yt,z,zt) (compute_acc(y,yt)*n+compute_acc(z,zt)*m)/(n+m);
        compute_overll_nmi = @(y,yt,z,zt) compute_mutual_info(blkdiag(y,z),blkdiag(yt,zt));
        
        y = generate_random_labels(n,K);
        z = generate_random_labels(m,L);
        ny = sum(y);
        nz = sum(z);
        P = prefactor*B*(log(m*n)^al)/sqrt(m*n);
        % P = prefactor*B*log(m*n)/sqrt(m*n);
        A = genSBM3(P,y,z);
         
        opts = {'ell1', true};%, 'type','US'};
     
        for mtd = 1:nReg
            tic,          
             [ yt, zt ] = biSpeClust( A, K, L, 'RegConst', RegVec(mtd), opts{:});
             yt = label_vec2mat(yt);
             zt = label_vec2mat(zt);
%             switch mtd
%                 case 1                
%                    [ yt, zt ] = biSpeClust( A, K, L, 'RegConst', RegVec(1), opts{:});
%                    yt = label_vec2mat(yt);
%                    zt = label_vec2mat(zt);
%                 case 2
%                    [ yt, zt ] = biSpeClust( A, K, L, 'RegConst', RegVec(2), opts{:});
%                    yt = label_vec2mat(yt);
%                    zt = label_vec2mat(zt);
%                 case 3
%                    [ yt, zt ] = biSpeClust( A, K, L, 'RegConst', RegVec(3), opts{:});
%                    yt = label_vec2mat(yt);
%                    zt = label_vec2mat(zt);
%                 case 4
%                    [ yt, zt ] = biSpeClust( A, K, L, 'RegConst', RegVec(4), opts{:});
%                    yt = label_vec2mat(yt);
%                    zt = label_vec2mat(zt);     
%             end       
            dt(j,mtd) = toc;
        
            nmi(j,mtd,1) = compute_mutual_info(yt,y);
            nmi(j,mtd,2) = compute_mutual_info(zt,z);
            nmi(j,mtd,3) = compute_overll_nmi(y, yt, z, zt);
            acc(j,mtd,1) = compute_acc(y, yt);
            acc(j,mtd,2) = compute_acc(z, zt);
            acc(j,mtd,3) = compute_overll_acc(y, yt, z, zt);
        end
        
        
      
    end

    nmi_all(:,:,:,t) = nmi;
    acc_all(:,:,:,t) = acc;
    dt_all(:,:,t) = dt;
%     java.io.File.createTempFile(sprintf('%3d',t), 'temp');
end

%%
%fprintf('\nrunning times = %s \n', sprintf('%4.2f  ', mean(dt_all,2)) )

%%
nmi_avg = mean(nmi_all,4);
acc_avg = mean(acc_all,4);
dt_avg = mean(dt_all,3);
result_fname = strrep(sprintf('results_C%2.2f_a%2.2f_b%2.2f_T%d_K%d_L%d',prefactor, al,beta, T,K,L),'.','p');
%save(sprintf('%s.mat',result_fname))

%%
%load('results.mat')
title_str = sprintf('C = %2.2f, \\alpha = %2.2f',prefactor,al);
figure(1), clf,
%subplot(131),hold on
%colors = parula(4);
colors = get(gca,'ColorOrder');
markers = {'-.','-s',':x','--'};
h = [];
for i = 1:nReg
    %h(i) = plot_ci_bands(n0_vec, squeeze(nmi_all(:,i,3,:)), colors(i,:), @plot);
    %h(i) = semilogx(n0_vec, nmi_avg(:,i,1),'LineWidth',2,'color',colors(i,:)), hold on
    h(i) = plot(n0_vec, nmi_avg(:,i,3),  markers{mod(i-1,4) + 1}, ...
        'LineWidth',2,'color',colors(i,:)); hold on
end

lgd = legend(h, sprintfc('%2.2f',RegVec), ... %{'SC1','SC2','SC3', 'SC4'}, ...
    'location','SouthEast');
legend('boxoff')
%axis([0,N+M,0,1]);
xlabel('# of nodes per cluster')
%ylabel('NMI (Overall)')
ylabel('overall NMI')
%title(title_str,'FontWeight','Normal')
axis tight
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 4.5];
fig.PaperPositionMode = 'manual';
%print('-depsc2',sprintf('%s_nmi.eps',result_fname))

%%
figure(2), clf, 
for i = 1:4  
    %plot_ci_bands(n0_vec, squeeze(acc_all(:,i,3,:)), colors(i,:), @plot)    
    %semilogx(n0_vec, log(1-acc_avg(:,i,1)),'LineWidth',2,'color',colors(i,:)), hold on
    h(i) = plot(n0_vec, nmi_avg(:,i,2),  markers{i}, ...
        'LineWidth',2,'color',colors(i,:));hold on
end
%ylabel(' log miss. (overall)')
ylabel('NMI Column')
xlabel('# of nodes per cluster')
axis tight
%title(title_str,'FontWeight','Normal')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 4.5];
fig.PaperPositionMode = 'manual';
%print('-depsc2',sprintf('%s_miss.eps',result_fname))
%ylim([0,1])

% figure(3), clf,
% for i = 1:4
%     %semilogy(n0_vec, dt_avg(:,i)),  hold on
%     plot(n0_vec, dt_avg(:,i)),  hold on
% end
% ylabel('acc (overall)')
% axis tight
% title(sprintf('prefactor = %2.2f',prefactor))
