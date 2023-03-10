function [RMSE_S11,RMSE_S12,RMSE_S22] = fourier_calc_RMSE_PBC(Pnm_2_11,Qnm_2_11,Rnm_2_11,Snm_2_11,...
                                        Pnm_2_12,Qnm_2_12,Rnm_2_12,Snm_2_12,...
                                        Pnm_2_22,Qnm_2_22,Rnm_2_22,Snm_2_22,...
                                        Lx,Ly,PBC_i,d,n_mat,m_mat,...
                                        spec_plot_subsetx,spec_plot_subsety,subx,suby,topleftxycoord,...
                                        MD_S11,MD_S12,MD_S22)
    
    %initialize
    sigxx_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigxy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigyy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    N = length(MD_S11(:)); 

    for j = 1:spec_plot_subsetx 
        for k = 1:spec_plot_subsety
            
            % get sigma (vir) value at point
            % NOTE: For now we assume, subset is small enough to be
            % represented as a point and not weighted area average
    
            %find midpoint of subset
            %NOTE: coord sys is now dislo which is center of subset
            midpoint = topleftxycoord + [j*subx-0.5*subx,-k*suby+0.5*suby]; %formula verified Research Log12/20 pg 249
            xvalue = midpoint(1);
            yvalue = midpoint(2);
    
            %nonpiecewise--------------------------------------------------------------
            %
            dx = xvalue;
            dy = yvalue;
            
            sigxx_2(k,j) = sigxx_2(k,j) +  sum(Pnm_2_11.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Qnm_2_11.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       +  Rnm_2_11.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Snm_2_11.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       ,'all');
            sigxy_2(k,j) = sigxy_2(k,j) +  sum(Pnm_2_12.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Qnm_2_12.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       +  Rnm_2_12.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Snm_2_12.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       ,'all');
            sigyy_2(k,j) = sigyy_2(k,j) + sum(Pnm_2_22.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                   +  Qnm_2_22.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                   +  Rnm_2_22.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                   +  Snm_2_22.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                   ,'all');
            
            %find summation terms
            for add_i = 1:PBC_i
    
			    dxR = xvalue-( add_i*d);
                dxL = xvalue-(-add_i*d);
			    dy = yvalue;
    
                %dislocations contributions for left dislocations
                dx = xvalue+add_i*d;
                dy = yvalue;
                sigxx_2(k,j) = sigxx_2(k,j) + sum(Pnm_2_11.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Qnm_2_11.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           +  Rnm_2_11.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Snm_2_11.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           ,'all');
                sigxy_2(k,j) = sigxy_2(k,j) + sum(Pnm_2_12.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Qnm_2_12.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           +  Rnm_2_12.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Snm_2_12.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           ,'all');
                sigyy_2(k,j) = sigyy_2(k,j) + sum(Pnm_2_22.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Qnm_2_22.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           +  Rnm_2_22.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Snm_2_22.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           ,'all');
                
    
    
                %dislocations contributions for right dislocations
                dx = xvalue-add_i*d;
                dy = yvalue;
    
		        sigxx_2(k,j) = sigxx_2(k,j) + sum(Pnm_2_11.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Qnm_2_11.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       +  Rnm_2_11.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Snm_2_11.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       ,'all');
		        sigxy_2(k,j) = sigxy_2(k,j) + sum(Pnm_2_12.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Qnm_2_12.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       +  Rnm_2_12.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Snm_2_12.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       ,'all');
		        sigyy_2(k,j) = sigyy_2(k,j) + sum(Pnm_2_22.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Qnm_2_22.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       +  Rnm_2_22.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Snm_2_22.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       ,'all');
            end %add_i = 1:PBC_i
            %end non piecewise-------------------------------------------------------------
            
        end %j = 1:subsety
    end %k = 1:subsetx

    %calculate RMSE between MD and reconstructed
    %bc RMSE is for minimization but the code is more maximization,
    % we will maximize the negative fitness which is effectively
    % minimization
    %RMSE will always produce a 0 or positive value, so the max of negative
    %should work
    RMSE_S11 = -sqrt(sum((sigxx_2-MD_S11).^2,'all')/N);
    RMSE_S12 = -sqrt(sum((sigxy_2-MD_S12).^2,'all')/N);
    RMSE_S22 = -sqrt(sum((sigyy_2-MD_S22).^2,'all')/N);
%     toc

%     H(1,i) = -RMSE_S11;
%     H(2,i) = -RMSE_S12;
%     H(3,i) = -RMSE_S22;

end %function