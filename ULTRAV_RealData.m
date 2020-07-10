
clear all;
close all;

% rng(5,'twister')
working_dir = '/home/tales/tmp/ULTRA_clean';

addpath(genpath(working_dir))
% addpath(genpath('~/Dropbox/matlab/HyperspectralCode/'))
SNR = 30;

% 
% cd /home/tales/MEGA/NovasImagens/Jasper_Ridge
cd(strcat(working_dir, '/Imgs/Jasper_Ridge'))

load jasperRidge2_R198.mat
% Y - HI in matrix form
bands = [1:99,105:137,146:198];
Y = Y(bands,:);
Y = Y/max(max(Y));
invRC = 1; % invert row-column
Y_cube = matrixToHCube(Y,nRow, nCol, invRC);
R = 4;
% M0 = vca(Y,'Endmembers',R);
load  t_jasper_em.mat 
M0 = M0(bands,:);
% M0 = find_endm(Y,R,'nfindr');

v = [57,30,20]; % RGB bands
rgb(:,:,1) = imadjust(rescale(Y_cube(:,:,v(2)),1));
rgb(:,:,2) = imadjust(rescale(Y_cube(:,:,v(1)),1));
rgb(:,:,3) = imadjust(rescale(Y_cube(:,:,v(3)),1));
figure, imshow(rgb); % display RGB image


% Samson data:

% cd /home/tales/MEGA/NovasImagens/Samson/
% load samson_1.mat
% % Y - HI in matrix form
% Y = V;
% Y = Y/max(max(Y));
% invRC = 1; % invert row-column
% Y_cube = matrixToHCube(Y,nRow, nCol, invRC);
% R = 3;
% %M0 = vca(Y,'Endmembers',R);
% load t_endmembers.mat


% cd /home/tales/MEGA/NovasImagens/Urban/
% load Urban_F210.mat
% Y = Y(SlectBands,:);
% %R=4;
% 
% invRC = 1; % invert row-column
% Y_cube = matrixToHCube(Y,nRow, nCol, invRC);
% [m,n,L] = size(Y_cube);
% % data_r = reshape(Y_cube,m*n,L);
% % r = data_r';
% % R=4;
% % M0 = vca(r,'Endmembers',R);
% load groundTruth_Urban_end6/end6_groundTruth.mat
% M0 = M;
% R = size(M,2);
% % Y_cube = Y_cube(50:250,80:end-20,:)/max(max(max(Y_cube)));
% % Y_cube = Y_cube(1:50,1:50,:)/max(max(max(Y_cube)));
% Y_cube = Y_cube/max(max(max(Y_cube)));
% 
% Y = hCubeToMatrix(Y_cube, invRC);
% [nRow,nCol, L] = size(Y_cube);


% cd /home/tales/MEGA/NovasImagens/
% load wash_dc_subim1.mat
% Y_cube = im;
% R=4;
% 
% invRC = 1; % invert row-column
% Y_cube = 0.9*Y_cube(1:80,1:100,:)/max(max(max(Y_cube)));
% Y = hCubeToMatrix(Y_cube, invRC);
% %  M0 = vca(Y,'Endmembers',R);
%  M0 = find_endm(Y,R,'nfindr');
% [nRow,nCol, L] = size(Y_cube);
% 


v = [57,30,20]; % RGB bands

[m,n,L] = size(Y_cube);

rgb(:,:,1) = imadjust(rescale(Y_cube(:,:,v(1)),1));
rgb(:,:,2) = imadjust(rescale(Y_cube(:,:,v(2)),1));
rgb(:,:,3) = imadjust(rescale(Y_cube(:,:,v(3)),1));
figure, imshow(rgb); % display RGB image

data_r = reshape(Y_cube,m*n,L);

% load endmembers_houston % load initial reference endmembers
% R = size(S0,2);

N = m*n;



% 
% materials{1} = 'vegetation'; % identify materials
% materials{2} = 'red roofs';
% materials{3} = 'concrete';
% materials{4} = 'asphalt';

r = data_r';

figure; plot(M0)
% Fully Constrained Least Squares Unmixing (FCLSU)
% M0 = vca(Y,'Endmembers',R);
disp('FCLSU...')
tic
A_FCLSU = FCLSU(data_r',M0)';
toc
A_FCLSU = A_FCLSU';

R_FCLS = M0*A_FCLSU';

% rmse_r_FCLS = RMSEAndSTDForMatrix(r(:), R_FCLS(:));
rmse_r_FCLS  = (norm(r(:)- R_FCLS(:))^2)/(N*L);

AAA = matrixToHCube(A_FCLSU',nRow,nCol,1);
figure;
for i=1:R
    subplot(2,R,i);
    imagesc(AAA(:,:,i),[0,1])
end
colormap(jet)

% Scaled version of the (partially) Constrained Least Squares (SCLSU)

disp('S-CLSU...')

tic
[A_SCLSU, psis_SCLSU] = SCLSU(Y_cube, M0);
toc

S_SCLSU = zeros(L,R,N);

parfor i = 1:m*n
   S_SCLSU(:,:,i) = M0*diag(psis_SCLSU(:,i)); 
end

% S_SCLSU = row2col_lexico_order(S_SCLSU,m,n);
acos_r_SCLS = 0;
acos_r_FCLS = 0;
R_SCLS = zeros(size(r));

for ll=1:N
    R_SCLS(:,ll) = squeeze(S_SCLSU(:,:,ll))*A_SCLSU(:,ll);
    acos_r_SCLS = acos_r_SCLS + acos(r(:,ll)'*R_SCLS(:,ll)/(norm(r(:,ll))*norm(R_SCLS(:,ll))));
    acos_r_FCLS = acos_r_FCLS + acos(r(:,ll)'*R_FCLS(:,ll)/(norm(r(:,ll))*norm(R_FCLS(:,ll))));
end
acos_r_SCLS = acos_r_SCLS/N;
acos_r_FCLS = acos_r_FCLS/N;
% rmse_r_SCLS = RMSEAndSTDForMatrix(r(:), R_SCLS(:));
rmse_r_SCLS  = (norm(r(:)- R_SCLS(:))^2)/(N*L);

AAA = matrixToHCube(A_SCLSU,nRow,nCol,1);
% figure;
for i=1:R
    subplot(2,R,i+R);
    imagesc(AAA(:,:,i),[0,1])
end

colormap(jet)


%% ULTRA-V

% Initializing A0 and M0 tensors:
A0_tensor = matrixToHCube(A_SCLSU,m,n,1);
% A0_tensor = matrixToHCube(A_FCLSU',m,n,1);
count = 1;
% for n1=1:m
%     for n2=1:n

for n2=1:n
    for n1=1:m
       M0_tensor(n1,n2,:,:)= M0;% + rand(size(M0))*0.001;
       MM_tensor(n1,n2,:,:)= S_SCLSU(:,:,count);
        count = count + 1;
    end
end

% Regularization parameter search ranges:
%lambda_A_search = [0.3];
% lambda_A_search = [0.1];
lambda_A_search = [0.00001];
%lambda_M_search = [0.8];
lambda_M_search = [0.8];
varepsilon =1e-3;
maxItt = 40;



% %Jasper Ridge
varepsilon =1e-4;
lambda_A_search = 0.00001;
lambda_M_search = 0.01;

% Estimating tensor ranks:
rank_A = intialTensorRankEst(A0_tensor)
rank_M = intialTensorRankEst(MM_tensor)


Nruns = length(lambda_A_search)*length(lambda_M_search); 

FullTensor_Results = zeros(Nruns,6);
fprintf('Nruns = %d\n\n', Nruns);
verbose = 1;
R_tensor = zeros(size(Y_cube));
acos_r_tensor = 0;
count = 1;
for lambda_A = lambda_A_search
    for lambda_M = lambda_M_search
        tic
        [ Mt, At , Pt, Qt] = unmixing_lr_reg(Y_cube, lambda_M, lambda_A, rank_M, rank_A, A0_tensor,MM_tensor,maxItt,verbose, varepsilon);
        toc
        acos_r_tensor = 0;
        for ll=1:m
            for jj=1:n
                R_tensor(ll,jj,:) = squeeze(Mt(ll,jj,:,:))*squeeze(At(ll,jj,:));
                acos_r_tensor = acos_r_tensor + acos(squeeze(Y_cube(ll,jj,:))'*squeeze(R_tensor(ll,jj,:))/( norm(squeeze(Y_cube(ll,jj,:)))*norm(squeeze(R_tensor(ll,jj,:))) ) );
            end                    
        end
        acos_r_tensor = acos_r_tensor/N;
        %rmse_r = RMSEAndSTDForMatrix(data(:), R_tensor(:));
        rmse_r = (norm(Y_cube(:)- R_tensor(:))^2)/(N*L);

        FullTensor_Results(count, :) = [lambda_A,lambda_M, rank_A, rank_M, rmse_r, acos_r_tensor];
        FullTensor_Results(count, :)
        count = count + 1;

    end
end


% ff = [file_str ,'Tensor_Houston_results.mat'];
% save(ff, 'FullTensor_Results');

A_FCLSU_im = matrixToHCube(A_FCLSU', nRow, nCol, invRC);
A_SCLSU_im = matrixToHCube(A_SCLSU, nRow, nCol, invRC);
figure;
for i=1:R
subplot(3,R,i);
imagesc(A_FCLSU_im(:,:,i),[0,1]);
subplot(3,R,i+R);
imagesc(A_SCLSU_im(:,:,i),[0,1]);
subplot(3,R,i+2*R);
imagesc(At(:,:,i),[0,1]);
end
colormap(jet)



%% Results and plots:
A_FCLSU_im = matrixToHCube(A_FCLSU', nRow, nCol, invRC);
A_SCLSU_im = matrixToHCube(A_SCLSU, nRow, nCol, invRC);


fprintf('FCLS: RMSE_R = %f, SAM_R = %f \n', 1e2*rmse_r_FCLS, 1e2*acos_r_FCLS)
fprintf('SCLS: RMSE_R = %f, SAM_R = %f \n', 1e2*rmse_r_SCLS, 1e2*acos_r_SCLS)
fprintf('ULTRA-V: RMSE_R = %f, SAM_R = %f \n', 1e2*FullTensor_Results(1,end-1),1e2*FullTensor_Results(1,end))


FSize = 16;

fh = figure;
[ha, pos] = tight_subplot(4, 3, 0.01, 0.1, 0.1);
for i=1:4
    axes(ha(1 + (i-1)*3));
    imagesc(A_FCLSU_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(2 + (i-1)*3));
    imagesc(A_SCLSU_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(3 + (i-1)*3));
    imagesc(At(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])    
end
set(fh, 'Position', [0 0 400 400])
 axes(ha(1));



%Jasper data:

title('FCLS','interpreter','latex','FontSize',FSize)
ylabel('Trees','interpreter','latex','FontSize',FSize)
axes(ha(2));
title('SCLSU','interpreter','latex','FontSize',FSize)
axes(ha(3));
title('ULTRA-V','interpreter','latex','FontSize',FSize)
axes(ha(4));
ylabel('Soil','interpreter','latex','FontSize',FSize)
axes(ha(7));
ylabel('Water','interpreter','latex','FontSize',FSize)
axes(ha(10));
ylabel('Road','interpreter','latex','FontSize',FSize)


colormap(jet)

% 


%%
ff = figure;
[haa, pos] = tight_subplot(1, R, 0.01, 0.1, 0.1);
for i=1:R
    ttt = hCubeToMatrix(Mt(:,:,:,i),invRC);
    ttt = ttt - repmat(mean(ttt,2),1,N);
    [V,D] = eig(ttt*ttt');
%     ttt2 = matrixToHCube(V(:,end)'*ttt, nRow,nCol,invRC);
%      ttt2 = matrixToHCube((D(end,end)^(-0.5))*V(:,end)'*ttt, nRow,nCol,invRC);
    ttt2 = matrixToHCube((D(end-R+1:end,end-R+1:end)^(-0.5))*V(:,end-R+1:end)'*ttt, nRow,nCol,invRC);
%     ttt2 = ttt2;%/max(max(max(ttt2)));
    ttt2 = mean(ttt2,3);
    ttt2 = ttt2 + abs(min(min(min(ttt2))));
    axes(haa(i));
%     imagesc(ttt2,[0 0.065])
imagesc(ttt2)
    set(gca,'ytick',[],'xtick',[])
end
set(ff, 'Position', [0 0 900 200])
colormap(jet)

% Samson
% axes(haa(1));
% title('Water','interpreter','latex','FontSize',FSize)
% axes(haa(2));
% title('Tree','interpreter','latex','FontSize',FSize)
% axes(haa(3));
% title('Soil','interpreter','latex','FontSize',FSize)


% Jasper Ridge
axes(haa(1));
title('Trees','interpreter','latex','FontSize',FSize)
axes(haa(2));
title('Soil','interpreter','latex','FontSize',FSize)
axes(haa(3));
title('Water','interpreter','latex','FontSize',FSize)
axes(haa(4));
title('Road','interpreter','latex','FontSize',FSize)





%%
FSize = 14;
fff = figure; 
[haaa, pos] = tight_subplot(1, R, 0.04, 0.1, 0.1);
for i=1:R
    temp = hCubeToMatrix(Mt(:,:,:,i));
    axes(haaa(i));
%     plot(temp(:,1:15:end),'color', ones(1,3)*0.9)
    % plot(temp(:,randsample(N,100)),'color', ones(1,3)*0.9
    
    if i==1
%         plot(temp(:,randsample(N,N/2))*0.8,'color', ones(1,3)*0.9)
        plot(temp(:,randsample(N,200))*0.8,'color', ones(1,3)*0.9)
         hold on
    	plot(M0(:,i),'k')    
    elseif i==2
%         plot(temp(:,randsample(N,N/2))*0.5,'color', ones(1,3)*0.9)
        plot(temp(:,randsample(N,200))*0.6,'color', ones(1,3)*0.9)
         hold on
    	plot(M0(:,i)*0.8,'k')    
    else
        plot(temp(:,randsample(N,200)),'color', ones(1,3)*0.9)
         hold on
        plot(M0(:,i),'k')    
    end
%      plot(M0(:,i),'k') 
    ylim([-0.01 1.05])
    xlim([0 L+1])
%     set(gca,'ytick',[],'xtick',[])
end
set(fff, 'Position', [0 0 900 240])
% Samson
% axes(haaa(1));
% title('Water','interpreter','latex','FontSize',FSize)
% axes(haaa(2));
% title('Tree','interpreter','latex','FontSize',FSize)
% axes(haaa(3));
% title('Soil','interpreter','latex','FontSize',FSize)

axes(haaa(1));
title('Trees','interpreter','latex','FontSize',FSize)
axes(haaa(2));
title('Soil','interpreter','latex','FontSize',FSize)
axes(haaa(3));
title('Water','interpreter','latex','FontSize',FSize)
axes(haaa(4));
title('Road','interpreter','latex','FontSize',FSize)



%%
fff = figure;
[haa, pos] = tight_subplot(2, R, 0.01, 0.1, 0.1);
for i=1:R
    axes(haa(i));    
    imagesc(At(:,:,i))    
    set(gca,'ytick',[],'xtick',[])
    axes(haa(i+R));
    imagesc(Qt(:,:,i))
    set(gca,'ytick',[],'xtick',[])
end
% set(ff, 'Position', [0 0 900 200])
colormap(jet)

axes(haa(1));
ylabel('${\mathcal{A}}$','interpreter','latex','FontSize',16);
axes(haa(1+R));
ylabel('${\mathcal{Q}}$','interpreter','latex','FontSize',16);
set(fff, 'Position', [0 0 900 240])

axes(haa(1));
title('Trees','interpreter','latex','FontSize',FSize)
axes(haa(2));
title('Soil','interpreter','latex','FontSize',FSize)
axes(haa(3));
title('Water','interpreter','latex','FontSize',FSize)
axes(haa(4));
title('Road','interpreter','latex','FontSize',FSize)



%%
FSize = 14;
fff = figure;
[haa, pos] = tight_subplot(2, R, 0.04, 0.1, 0.1);
for i = 1:R
    tempM = hCubeToMatrix(Mt(:,:,:,i), invRC);
    tempP = hCubeToMatrix(Pt(:,:,:,i), invRC);
    if i==1
        tempM = tempM* 0.8;
        tempP = tempP* 0.8;
    elseif i==2
        tempM = tempM* 0.6;
        tempP = tempP* 0.6;
    end
    axes(haa(i));    
    plot(tempM(:,[1:100:end]))%,'color', ones(1,3)*0.9); 
    samp = [1:50:N];
    ymax = 1.1*max(max(tempM(:,samp)));    
    ylim([-0.05 ymax])
%     set(gca,'ytick',[],'xtick',[])
    set(gca,'xtick',[])
    axes(haa(i+R));    
    plot(tempP(:,[1:50:end]))%,'color', ones(1,3)*0.9)
    ylim([-0.1 ymax])
    set(gca,'xtick',[])
end
axes(haa(1));
ylabel('${\mathcal{M}}$','interpreter','latex','FontSize',16);
axes(haa(1+R));
ylabel('${\mathcal{P}}$','interpreter','latex','FontSize',16);
set(fff, 'Position', [0 0 900 240])

axes(haa(1));
title('Trees','interpreter','latex','FontSize',FSize)
axes(haa(2));
title('Soil','interpreter','latex','FontSize',FSize)
axes(haa(3));
title('Water','interpreter','latex','FontSize',FSize)
axes(haa(4));
title('Road','interpreter','latex','FontSize',FSize)

