%LD Lx is subject to change which means you must adjust spec_plot_xlength
%and ylength
%the 2.4 in subsetx and subsety could also be subject to change
function H=evaluation_parfor(Pnm_2_11_full,Qnm_2_11_full,Rnm_2_11_full,Snm_2_11_full,...
                        Pnm_2_12_full,Qnm_2_12_full,Rnm_2_12_full,Snm_2_12_full,...
                        Pnm_2_22_full,Qnm_2_22_full,Rnm_2_22_full,Snm_2_22_full,...
                        Lx,Ly,PBC_i,d,n_mat,m_mat,...
                        spec_plot_subsetx,spec_plot_subsety,subx,suby,topleftxycoord,...
                        MD_S11,MD_S12,MD_S22)



% PQRS = population
% n = chromosomes to be mutated
% pages = number of chromosomes
[term_n,term_m,pages]=size(Pnm_2_11_full);

H=zeros(3,pages); % S11 S12 S22 x number of chromosomes
H_S11=zeros(1,pages); % S11 x number of chromosomes
H_S12=zeros(1,pages); % S12 x number of chromosomes
H_S22=zeros(1,pages); % S22 x number of chromosomes

%Old code
% %import weighted area MD 
% %MD Lx = 50.4;
% MD_S11 = readmatrix(strcat(MDpath,'_S11'));
% MD_S12 = readmatrix(strcat(MDpath,'_S12'));
% MD_S22 = readmatrix(strcat(MDpath,'_S22'));
% 
% % %for test only
% % [testx,testy] = size(MD_S11);
% % MD_S11 = MD_S11(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
% % MD_S12 = MD_S12(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
% % MD_S22 = MD_S22(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
% 
% %create n x m matrix
% m_mat = zeros(term_m,term_n); %y=row x=col
% n_mat = zeros(term_m,term_n); %y=row x=col
% m_row = 1:term_m;
% n_row = 1:term_n;
% for len = 1:term_n %for every column
%     m_mat(:,len) = m_row';
% end %len = 1:n
% for len = 1:term_m %for every row
%     n_mat(len,:) = n_row;
% end %len = 1:m
% 
% 
% %Lx and Ly of reconstructed GB is set manually
% 
% %the following is defined above--------------------------------------------
% %define new spectral plot values-------------------------------------------
% spec_plot_xlength = 50.4;
% spec_plot_ylength = 50.4;
% spec_plot_subsetx = 2*spec_plot_xlength/2.4;
% spec_plot_subsety = 2*spec_plot_xlength/2.4;
% subx = 2*spec_plot_xlength/spec_plot_subsetx; % xdist value of each subsetx
% suby = 2*spec_plot_ylength/spec_plot_subsety; % ydist value of each subsety
% subArea = subx*suby;
% 
% topleftxycoord = [-spec_plot_xlength ,spec_plot_ylength];

%reconstruct stress field of GB using constants and fourier series
parfor i = 1:pages

    [H_S11(i),H_S12(i),H_S22(i)] = fourier_calc_RMSE_PBC(Pnm_2_11_full(:,:,i),Qnm_2_11_full(:,:,i),Rnm_2_11_full(:,:,i),Snm_2_11_full(:,:,i),...
                                            Pnm_2_12_full(:,:,i),Qnm_2_12_full(:,:,i),Rnm_2_12_full(:,:,i),Snm_2_12_full(:,:,i),...
                                            Pnm_2_22_full(:,:,i),Qnm_2_22_full(:,:,i),Rnm_2_22_full(:,:,i),Snm_2_22_full(:,:,i),...
                                            Lx,Ly,PBC_i,d,n_mat,m_mat,...
                                            spec_plot_subsetx,spec_plot_subsety,subx,suby,topleftxycoord,...
                                            MD_S11,MD_S12,MD_S22);

end %for i = 1:pages

H(1,:) = H_S11;
H(2,:) = H_S12;
H(3,:) = H_S22;