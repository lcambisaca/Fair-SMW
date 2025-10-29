close all; clear; clc;
format long
anz_runs = 1;
n_range = 1000:1000:10000;
alg5_gap = zeros(1,10)
alg3_gap = zeros(1,10)
alg6RW_gap = zeros(1,10)
alg6SYM_gap = zeros(1,10)



k = 2;
h = 2;

a = 0.000060;
b = 0.000050;
c = 0.000005;
d = 0.000002;
% a = 0.9999;
% b = 0.9998;
% c = 0.9997;
% d = 0.9996; 



time_alg1 = zeros(anz_runs,length(n_range)); %Nomalized-SC
time_alg3 = zeros(anz_runs,length(n_range)); %S-Fair-SC
time_alg5 = zeros(anz_runs,length(n_range)); %W+n
time_alg6Rw = zeros(anz_runs,length(n_range)); %RW-SVM-SC: G = D\W + n * eye(n)
time_alg6Sym = zeros(anz_runs,length(n_range)); %SYM-SVM-SC: G = (sqrtD\W)/sqrtD +  n * eye(n)

error_alg1 = zeros(anz_runs,length(n_range));
error_alg3 = zeros(anz_runs,length(n_range));
error_alg5 = zeros(anz_runs,length(n_range));
error_alg6Rw = zeros(anz_runs,length(n_range));
error_alg6Sym = zeros(anz_runs,length(n_range));

eigs_alg1 = zeros(anz_runs,length(n_range));
eigs_alg3 = zeros(anz_runs,length(n_range));
eigs_alg5 = zeros(anz_runs,length(n_range));
eigs_alg6Rw = zeros(anz_runs,length(n_range));
eigs_alg6Sym = zeros(anz_runs,length(n_range));


for mmm = 1:length(n_range)
    n = n_range(mmm);

    fprintf('-----------------n = %d --------------------------------\n', n);

    block_sizes = (n/(k*h))*ones(1,k*h);
    sensitive0 = zeros(n,1);
    labels = zeros(n,1);   % <-- keep an untouched copy per n

    for yyy = 1:k
        for zzz = 1:h
            idx1 = ((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1);
            idx2 = ((n/k)*(yyy-1)+(n/(k*h))*zzz);
            sensitive0(idx1:idx2) = zzz;
            labels0(idx1:idx2) = yyy;
        end
    end

    for ell = 1:anz_runs
        % rng(1000 + n + ell); 


        fprintf('-----------------run = %d --------------------------------\n', ell);

        [W, D, F] = generate_connected_SBM(n,a,b,c,d,k,h,block_sizes,sensitive0);

        n1 = size(W,1)
        edges1 = nnz(W);
        density = edges1/(n1*(n1-1) );
        fprintf('    density = %.6f\n', density);
       
        tstart1 = tic;
        [labelsalg1 , t1] = alg1(W, D, k);
        time_alg1(ell,mmm) = toc(tstart1);
        eigs_alg1(ell,mmm) = t1;
        error_alg1(ell,mmm)=clustering_accuracy(labels,labelsalg1);    

        tstart3 = tic;
        [labelsalg3 , t3,gap3] = alg3(W, D, F, k);
        time_alg3(ell,mmm) = toc(tstart3);
        eigs_alg3(ell,mmm) = t3;
        error_alg3(ell,mmm)=clustering_accuracy(labels,labelsalg3);   
        alg3_gap(ell,mmm) = gap3


  
        tstart5 = tic;
        [labelsalg5 , t6,gap5] = alg5(W, D, F, k);
        time_alg5(ell,mmm) = toc(tstart5);
        eigs_alg5(ell,mmm) = t6;
        error_alg5(ell,mmm)=clustering_accuracy(labels,labelsalg5);
        alg5_gap(ell,mmm) = gap5


        tstart6Rw = tic;
        [labelsalg6Rw , t7,gap6rw] = alg6Rw(W, D, F, k);
        time_alg6Rw(ell,mmm) = toc(tstart6Rw);
        eigs_alg6Rw(ell,mmm) = t7;
        error_alg6Rw(ell,mmm)=clustering_accuracy(labels,labelsalg6Rw);
        alg6RW_gap(ell,mmm) = gap6rw



        tstart6Sym = tic;
        [labelsalg6Sym , t8,gap6sym] = alg6Sym(W, D, F, k);
        time_alg6Sym(ell,mmm) = toc(tstart6Sym);
        eigs_alg6Sym(ell,mmm) = t8;
        error_alg6Sym(ell,mmm)=clustering_accuracy(labels,labelsalg6Sym);
        alg6SYM_gap(ell,mmm) = gap6sym


    end
end

difftime1 = time_alg1 - eigs_alg1;
difftime3 = time_alg3 - eigs_alg3;
difftime5 = time_alg5 - eigs_alg5;
difftime6Rw = time_alg6Rw - eigs_alg6Rw;
difftime6Sym = time_alg6Sym - eigs_alg6Sym;

% set default sizes for figures:
ulesfontsize = 15;
set(0, 'DefaultAxesFontSize', ulesfontsize);
set(0, 'DefaultTextFontSize', ulesfontsize);
set(0, 'DefaultUIControlFontSize', ulesfontsize);
set(0,'DefaultLineMarkerSize',ulesfontsize);
set(0,'DefaultLineLineWidth',2.5)
set(gcf, 'PaperPosition', [0 0 10 7.5])
set(gcf, 'PaperSize', [10 7.5]);

% --- Runtime Plot ---
n2_curve = n_range.^2;     % <-- match n_range
scale_factor = mean(mean(time_alg5,1)) / mean(n2_curve);
scaled_n2_curve = scale_factor * n2_curve;

figure;clf;
hold on
plot(n_range, median(error_alg1,1), 'bx-');   % blue x
% plot(n_range, median(error_alg2,1), 'mo-');   % magenta circle
plot(n_range, median(error_alg3,1), 'rd-');   % yellow diamond
% plot(n_range, median(error_alg4Rw,1), 'r^-'); % red triangle up
% plot(n_range, median(error_alg4Sym,1), 'k+-');% black plus
plot(n_range, median(error_alg5,1), 'gs-');   % green square
plot(n_range, median(error_alg6Sym,1), 'c*-'); % cyan star
plot(n_range, median(error_alg6Rw,1), 'mo-');% magenta | 
hold off
legend({'SC','S-Fair-SC','AFF-Fair-SMW','SYM-Fair-SNW','RW-Fair-SMW'}, 'Location','northwest', 'FontSize',9)
xlabel('n')
ylabel('Error')
ylim([0,1])
title(strcat('Error Rate of Algorithms using SBM Model' , ' (h=',num2str(h),', k=',num2str(k) , ')'),'FontWeight','normal')
grid on

figure;clf;
hold on
plot(n_range, mean(time_alg1,1), 'bx-');       
plot(n_range, mean(time_alg3,1), 'rd-');       
plot(n_range, mean(time_alg5,1), 'gs-');       
plot(n_range, mean(time_alg6Sym,1), 'c*-');    
plot(n_range, mean(time_alg6Rw,1), 'mo-');    
plot(n_range, scaled_n2_curve, 'k--');
hold off
legend({'SC','S-Fair-SC','AFF-SVM-SC','SYM-SVM-SC','RW-SVM-SC','O(n^2)'}, 'Location','northwest', 'FontSize',9)
xlabel('n')
ylabel('Running time (s)')
xticks(n_range)            % <-- force exact n labels
title(strcat('Algorithms Run Time Using SBM Model ', ' (h=',num2str(h),', k=',num2str(k),')'), 'FontWeight','normal')
grid on

% --- Eigs Time ---
figure;clf;
hold on
plot(n_range,mean(eigs_alg1,1),'bx-')
plot(n_range,mean(eigs_alg3,1),'rd-')
plot(n_range,mean(eigs_alg5,1),'gs-')
plot(n_range,mean(eigs_alg6Sym,1),'c*-')
plot(n_range,mean(eigs_alg6Rw,1),'mo-')
hold off
legend({'SC','S-Fair-SC','AFF-SVM-SC','SYM-SVM-SC','RW-SVM-SC'}, 'Location','northwest', 'FontSize',9)
xlabel('n')
ylabel('Running time (s)')
xticks(n_range)            % <-- force exact n labels
title(strcat('Algorithms Eigs Run Time Using SBM Model ', ' (h=',num2str(h),', k=',num2str(k),')'), 'FontWeight','normal')
grid on

% --- Runtime excluding eigs ---
figure;clf;
hold on
plot(n_range,mean(difftime1,1),'bx-')
plot(n_range,mean(difftime3,1),'rd-')
plot(n_range,mean(difftime5,1),'gs-')
plot(n_range,mean(difftime6Sym,1),'c*-')
plot(n_range,mean(difftime6Rw,1),'mo-')
hold off
legend({'SC','S-Fair-SC','AFF-SVM-SC','SYM-SVM-SC','RW-SVM-SC'}, 'Location','northwest', 'FontSize',9)
xlabel('n')
ylabel('Running time (s)')
xticks(n_range)            % <-- force exact n labels
title(strcat('Algorithms Run-time excluding eigs Using SBM Model ', ' (h=',num2str(h),', k=',num2str(k),')'), 'FontWeight','normal')
grid on
