function plotfourier_reconstructed(Pnm_2_11,Qnm_2_11,Rnm_2_11,Snm_2_11,...
                                        Pnm_2_12,Qnm_2_12,Rnm_2_12,Snm_2_12,...
                                        Pnm_2_22,Qnm_2_22,Rnm_2_22,Snm_2_22,...
                                        Lx,Ly,PBC_i,d,n_mat,m_mat,...
                                        spec_plot_xlength,spec_plot_ylength,spec_plot_subsetx,spec_plot_subsety,subx,suby,topleftxycoord,...
                                        MD_S11,MD_S12,MD_S22,...
                                        gen,destdir)
    
    %initialize
    sigxx_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigxy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigyy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);

    parfor j = 1:spec_plot_subsetx 
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

    %plot and save differences----------------------------------------------
    Opt_S11 = sigxx_2;
    Opt_S12 = sigxy_2;
    Opt_S22 = sigyy_2;
    %%
    [testx,testy] = size(Opt_S11);
    Opt_S11 = Opt_S11(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
    Opt_S12 = Opt_S12(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
    Opt_S22 = Opt_S22(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
    
    %%
    diff_S11 = MD_S11-Opt_S11;
    diff_S12 = MD_S12-Opt_S11;
    diff_S22 = MD_S22-Opt_S11;
    
    %%
    destdir_gen = strcat(destdir,'gen',num2str(gen),'\');
    mkdir(destdir_gen);
    
    %optimal
    figName= 'Optimal S11 Contours';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),Opt_S11,'CDataMapping','scaled')
    caxis([min(Opt_S11(:)),max(Opt_S11(:))])
    colorbar
    axis equal
    saveas(f, [destdir_gen,figName,'.png']);
    
    figName= 'Optimal S12 Contours';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),Opt_S12,'CDataMapping','scaled')
    caxis([min(Opt_S12(:)),max(Opt_S12(:))])
    colorbar
    axis equal
    saveas(f, [destdir_gen,figName,'.png']);
    
    figName= 'Optimal S22 Contours';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),Opt_S22,'CDataMapping','scaled')
    caxis([min(Opt_S22(:)),max(Opt_S22(:))])
    colorbar
    axis equal
    saveas(f, [destdir_gen,figName,'.png']);
    
    %difference
    figName= 'diff S11 Contours';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),diff_S11,'CDataMapping','scaled')
    caxis([min(diff_S11(:)),max(diff_S11(:))])
    colorbar
    axis equal
    saveas(f, [destdir_gen,figName,'.png']);
    
    figName= 'diff S12 Contours';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),diff_S12,'CDataMapping','scaled')
    caxis([min(diff_S12(:)),max(diff_S12(:))])
    colorbar
    axis equal
    saveas(f, [destdir_gen,figName,'.png']);
    
    figName= 'diff S22 Contours';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),diff_S22,'CDataMapping','scaled')
    caxis([min(diff_S22(:)),max(diff_S22(:))])
    colorbar
    axis equal
    saveas(f, [destdir_gen,figName,'.png']);

end %end function