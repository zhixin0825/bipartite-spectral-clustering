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
n0 = 500;
m = L*n0;
n = K*n0;

prefactor = 1;
al = .5;
P = prefactor*B*(log(m*n)^al)/sqrt(m*n);
beta = 1;

tau_len = 20;
tau_vec = logspace(-1,1,tau_len);
%%
nmtds = 3;

T = 2;
err_all = zeros(tau_len, nmtds, T);


for t = 1:T  % This is parallel for; use with caution
    %nmi = zeros(n0_len,4,2);
    sprintf('.');
    err = zeros(tau_len, nmtds);
  
    y = generate_random_labels(n,K);
    z = generate_random_labels(m,L);
    A = genSBM3(P,y,z);
         
    %M = y * P * z';
    d1 = y*P*sum(z)';
    d2 = z*P'*sum(y)';
    dmax1 = max(d1);
    dmax2 = max(d2);
    normM = compMeanNorm(y,z,P); %norm(M);
    
    yP = y*P;
    zPt = z*P';
    err(1,1) = compDevFromMean(A,y,z,yP,zPt) / normM; % norm(A-M)/normM;
    for j = 1:tau_len
        fprintf('.')
        err(j,1) = err(1,1);
        
        Are = regularizeAdj(A,'tau',tau_vec(j),'ell1',false);
        err(j,2) =  compDevFromMean(Are,y,z,yP,zPt) / normM; %norm(Are-M)/normM;
        
        Are = regularizeAdj(A,'tau',tau_vec(j),'dmax1', dmax1, 'dmax2', dmax2, 'ell1', false);
        err(j,3) =  compDevFromMean(Are,y,z,yP,zPt) / normM; %norm(Are-M)/normM;
    end
    fprintf('\n')
    
    err_all(:,:,t) = err;
end

%%
result_fname = strrep(sprintf('results_C%2.2f_a%2.2f_b%2.2f_T%d_K%d_L%d',prefactor, al,beta, T,K,L),'.','p');
%save(sprintf('%s.mat',result_fname))

%%
err_avg = mean(err_all,3);
%load('results.mat')
title_str = sprintf('C = %2.2f, \\alpha = %2.2f',prefactor,al);
figure(1), clf,
colors = get(gca,'ColorOrder');
markers = {'-.','--s',':x','--'};
h = [];
for i = 1:nmtds
    %h(i) = plot_ci_bands(n0_vec, squeeze(nmi_all(:,i,3,:)), colors(i,:), @plot);
    %h(i) = semilogx(n0_vec, nmi_avg(:,i,1),'LineWidth',2,'color',colors(i,:)), hold on
    h(i) = semilogx(tau_vec, err_avg(:,i),  markers{i}, ...
        'LineWidth',2,'color',colors(i+3,:)); hold on
end

lgd = legend(h, {'No reg.', 'data-driv. Reg.', 'dmax Reg.'},'Position',[0.65 .55 0.2 0.2]);
legend('boxoff')
%axis([0,N+M,0,1]);
xlabel('Regularization threshold ($\tau$)','interpreter','latex')
%ylabel('NMI (Overall)')
ylabel('Relative operator norm error')
title(title_str,'FontWeight','Normal')
%axis tight
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 4.5];
fig.PaperPositionMode = 'manual';
%print('-dpng','-r600',sprintf('%s_concent.png',result_fname))
%print('-depsc',sprintf('%s_concent.eps',result_fname))

