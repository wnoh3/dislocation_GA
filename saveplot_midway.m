%only run this code if you stop the ga_parfor code during the gen loop
% DO NOT CLEAR ANY VARIABLES FROM PREVIOUS RUN
fprintf("GA Finished.\n");

%% save score plots
figure(1);
figName= 'maxAvgScores_S11';
saveas(figure(1), [destdir,figName,'.png']);

figure(2);
figName= 'maxAvgScores_S12';
saveas(figure(2), [destdir,figName,'.png']);

figure(3);
figName= 'maxAvgScores_S22';
saveas(figure(3), [destdir,figName,'.png']);

Max_fitness_value_S11=max(K_11(:,2)) %find best of best
Max_fitness_value_S12=max(K_12(:,2)) %find best of best
Max_fitness_value_S22=max(K_22(:,2)) %find best of best

Optimal_solution_P_11 = P_11(:,:,1); % Best chromosome
Optimal_solution_Q_11 = Q_11(:,:,1); % Best chromosome
Optimal_solution_R_11 = R_11(:,:,1); % Best chromosome
Optimal_solution_S_11 = S_11(:,:,1); % Best chromosome

Optimal_solution_P_12 = P_12(:,:,1); % Best chromosome
Optimal_solution_Q_12 = Q_12(:,:,1); % Best chromosome
Optimal_solution_R_12 = R_12(:,:,1); % Best chromosome
Optimal_solution_S_12 = S_12(:,:,1); % Best chromosome

Optimal_solution_P_22 = P_22(:,:,1); % Best chromosome
Optimal_solution_Q_22 = Q_22(:,:,1); % Best chromosome
Optimal_solution_R_22 = R_22(:,:,1); % Best chromosome
Optimal_solution_S_22 = S_22(:,:,1); % Best chromosome

%%
%plot using optimal solution
[Opt_S11,Opt_S12,Opt_S22] = plot_fourier_isolated(Optimal_solution_P_11,Optimal_solution_Q_11,Optimal_solution_R_11,Optimal_solution_S_11,...
                                                    Optimal_solution_P_12,Optimal_solution_Q_12,Optimal_solution_R_12,Optimal_solution_S_12,...
                                                    Optimal_solution_P_22,Optimal_solution_Q_22,Optimal_solution_R_22,Optimal_solution_S_22,...
                                                    Lx,Ly,n_mat,m_mat,...
                                                    spec_plot_xlength,spec_plot_ylength,spec_plot_subsetx,spec_plot_subsety,subx,suby,topleftxycoord,...
                                                    destdir);


