function varargout = UI(varargin)
clear;

    
global i_f_method i_lambda i_I i_R i_temporal i_seasonal i_periodic i_mode_Periodic i_mode_Seasonal i_I_Seasonal i_R_Seasonal i_I_Periodic i_R_Periodic i_struct i_order i_order_nnt i_s_angle i_s_corr i_defiend_periodic i_defiend_anomalies i_defiend_changes i_added i_streaming i_noise i_nonnegative i_sparse i_normal i_fixsign i_ns X; 
i_I=[];
i_R=[];
i_f_method={};
i_temporal=0;
i_seasonal=0;
i_periodic=0;
i_mode_Periodic=0;
i_mode_Seasonal=0;
i_R_Seasonal=0;
i_I_Seasonal=0;
i_I_Periodic=0;
i_R_Periodic=0;
i_struct='';
i_order=[];
i_order_nnt=0;
i_s_angle=0;
i_s_corr=0;
i_defiend_periodic=0;
i_defiend_seasonal=0;
i_defiend_anomalies=0;
i_defiend_changes=0;
i_streaming=0;
i_added=0;
i_noise=0;
i_nonnegative=0;
i_sparse=0;
i_sparse_param=0.5;
i_normal=0;
i_fixsign=0;
i_lambda=[1 1 1];
i_ns=0;
tabs={'Modes','PeriodWave','SeasonalEffects','Streaming','ChangePoints','Anomalies','Noise','NonNegative','Sparsity','Settings','Generate','About'};

if ~isdeployed 
    addpath('functions');
    addpath('tensor_toolbox_2.6');
end


% *********************************************************************************************************
% F U N C T I O N S
% *********************************************************************************************************
% Export tensor
function save_tensor(hObject, eventdata)
    [file,path,filterindex] = uiputfile({'*.mat','Tensor Toolbox Format (*.mat)';'*.mat','Dense MATLAB format (*.mat)';'*.csv','Comma Separated Values (*.csv)';'*.hdf','Hierarchical Data Format - HDF5 (*.hdf)'},'Save file name');
    if file
        filename=strcat(path,file);
        tobj=tab{findtab('Generate')};
        edm=findobj(tobj,'Tag','edit_multiline');
        line=get(edm, 'String');
        filetypes={'Tensor Toolbox','Dense MATLAB array','CSV','HDF5'};
        line{length(line)+1}=sprintf('Exporting %s format ... Please wait ...',filetypes{filterindex});
        set(edm, 'String',line);
        pause(0.5);
        if exist(filename, 'file') == 2
            delete(filename);
        end   
        switch filterindex
            case 1 % Tensor Toolbox Format (incl. sparse tensors)
                save(filename,'X');
            case  2 % Dense Matlab array
                if ~isa(X,'double')
                    Xd=double(X);
                    save(filename,'Xd');
                else
                    save(filename,'X');
                end    
            case  3 % CSV
               if isa(X,'sptensor')
                    csvwrite(filename,[X.subs X.vals]);
                else
                    Xs=sptensor(X);
                    csvwrite(filename,[Xs.subs Xs.vals]);
                end
            case  4 %HDF
                if isa(X,'sptensor')
                    dat=[X.subs X.vals];
                    h5create(filename,'/X',size(X));
                    h5write(filename,'/X',dat);   
                else
                    Xd=double(X);
                    h5create(filename,'/X',size(Xd));
                    h5write(filename,'/X',Xd);   
                end
        end
        line{length(line)+1}=sprintf('Saving to "%s" ... Done',filename);
        set(edm, 'String',line);
    end    
end    

% Generate Tensor
function generate_tensor(hObject, eventdata)
    tobj=tab{findtab('Generate')};
    edm=findobj(tobj,'Tag','edit_multiline');
    savebutton=findobj(tobj,'Tag','savebutton');

    % **************************************************************
    % initialization check
    % **************************************************************
    MaxP=0;
    MaxS=0;
    if i_order-i_order_nnt==1
        if i_seasonal && i_periodic
            MaxP=i_R(i_order_nnt+1); 
            MaxS=i_R(i_order_nnt+1);
            if i_defiend_periodic+i_defiend_seasonal~=MaxP
              errordlg(sprintf('You have defined "%i" factors in total for peridoic waves and seasonal effects while it should be equal to "%i" in total!',i_defiend_periodic+i_defiend_seasonal,i_R(i_order_nnt+1)));
              tgroup.SelectedTab = tab{findtab('PeriodWave')};
              return;
            end    
        end 
        if i_seasonal==1 && i_periodic==0
            MaxS=i_R(i_order_nnt+1);
            if i_defiend_seasonal~=MaxS
               errordlg(sprintf('You have defined "%i" factors in total for seasonal effects while it should be equal to "%i"!',i_defiend_seasonal,i_R(i_order_nnt+1)));
               tgroup.SelectedTab = tab{findtab('SeasonalEffects')};
               return;
            end   
        end
        if i_seasonal==0 && i_periodic==1
            MaxP=i_R(i_order_nnt+1);
            if i_defiend_periodic~=MaxP
               errordlg(sprintf('You have defined "%i" factors in total for periodic waves while it should be equal to "%i"!',i_defiend_periodic,i_R(i_order_nnt+1)));
              tgroup.SelectedTab = tab{findtab('PeriodWave')};
              return;              
            end 
        end
    elseif i_order-i_order_nnt==2
        MaxP=i_R(i_order_nnt+1);
        MaxS=i_R(i_order_nnt+2);
        if i_defiend_periodic~=MaxP
             errordlg(sprintf('You have defined "%i" factors in total for periodic waves while it should be equal to "%i"!',i_defiend_periodic,i_R(i_order_nnt+1)));
             tgroup.SelectedTab = tab{findtab('PeriodWave')};
             return;
        end
        if i_defiend_seasonal~=MaxS
             errordlg(sprintf('You have defined "%i" factors in total for seasonal effects while it should be equal to "%i"!',i_defiend_seasonal,i_R(i_order_nnt+1)));
             tgroup.SelectedTab = tab{findtab('SeasonalEffects')};
             return;
        end
    end    
   
    if i_defiend_changes
        tobj_changes=tab{findtab('ChangePoints')};
        data_changes_p=get(findobj(tobj_changes,'Tag','uitable_changepoints_p'),'Data');
        data_changes_s=get(findobj(tobj_changes,'Tag','uitable_changepoints_s'),'Data');
    
        for n=1:size(data_changes_p,1)
            if cell2mat(data_changes_p(n,1))>i_I(i_mode_Periodic) || cell2mat(data_changes_p(n,2))>i_I(i_mode_Periodic) 
                errordlg(sprintf('The start point for changepoint exceed the defined size of dimension. You should reset changepoints!'));
                tgroup.SelectedTab = tab{findtab('ChangePoints')};
                return;
            end     
        end    
        for n=1:size(data_changes_s,1)
            if cell2mat(data_changes_s(n,1))>i_I(i_mode_Seasonal) || cell2mat(data_changes_s(n,2))>i_I(i_mode_Seasonal) 
                 errordlg(sprintf('The start point for changepoint exceed the defined size of dimension. You should reset changepoints!'));
                tgroup.SelectedTab = tab{findtab('ChangePoints')};
                return;
            end     
        end    
        cell2mat(data_changes_s(:,1))
        cell2mat(data_changes_s(:,2))
    end
    
   % ************************
   % end
   % ************************
    line{1}=sprintf('********** SimTensor v1.0 - developed by Hadi Fanaee-T, INESC TEC *************');
    line{2}=sprintf('Generating %s-sctructured Tensor with order of %i and size of [%s]...',i_struct,i_order,strrep(num2str(i_I),'  ',' '));
    set(edm, 'String',line);
    if strcmp(i_struct,'CP')
        typ='Lambda';
        lambdastr=sprintf('[%s]',strrep(num2str(i_lambda),'  ',' '));
    else
        typ='G';
        lambdastr=sprintf('[%s] (%s)',strrep(num2str(i_lambda),'  ',' '),strrep(num2str(i_R),'  ',' '));
    end    
    line{length(line)+1}=sprintf(' %s = %s',typ,lambdastr);
    set(edm, 'String',line);
    line{length(line)+1}=sprintf(' Generating Factor matrices ...');
    set(edm, 'String',line)
    for n=1:i_order_nnt
  
       line{length(line)+1}=sprintf('  Mode #%i (Non-temporal): Generate "%i" factors with "%s" distribution, length of "%i" and noise "%f"',n,i_R(n),i_f_method{n},i_I(n),i_ns);
       set(edm, 'String',line);
        if i_s_angle
            line{length(line)+1}=sprintf('  Mode #%i Applying Angle of "%f" between factors',n,i_s_angle);
            set(edm, 'String',line);
            U0{n}=SimAngle(i_I(n),i_R(n),i_s_angle,eval(strcat('@',i_f_method{n})));
        else
            U0{n}=SimFacGen(i_f_method{n},i_I(n),i_R(n),i_ns);
        end
        if i_s_corr
             line{length(line)+1}=sprintf('  Mode #%i Applying Correlation of "%f" between factors',n,i_s_corr);
             set(edm, 'String',line);
             U0{n}=SimCol(U0{n},i_s_corr);
        end 
    end
    
    if i_periodic
        
        line{length(line)+1}=sprintf('  Mode #%i (Perioidc waves)',i_mode_Periodic);
        set(edm, 'String',line);
        tobj_periodic=tab{findtab('PeriodWave')};
        data_periodic=get(findobj(tobj_periodic,'Tag','uitable_periodic_table'),'Data');

        %i_p_method = strjoin(data_periodic(:,1),'|');
        %i_p_nwaves = strrep(strjoin(cellstr(num2str(cell2mat(data_periodic(:,2)))),'|'),' ','');
        %i_p_freq = strrep(strjoin(cellstr(num2str(cell2mat(data_periodic(:,3)))),'|'),' ','');
        for n=1:size(data_periodic,1)
            line{length(line)+1}=sprintf('   Periodic wave factor #%i with "%s" function, "%i" waves with frequency "%i" , and noise "%f"',n,strjoin(data_periodic(n,1)),cell2mat(data_periodic(n,2)),cell2mat(data_periodic(n,3)),i_ns );
            set(edm, 'String',line);
            U0{i_mode_Periodic}(:,n)=SimPeriodWave(strjoin(data_periodic(n,1)),i_I_Periodic,cell2mat(data_periodic(n,2)),cell2mat(data_periodic(n,3)),i_ns);
        end
    end   
    
    if i_seasonal
        
        if i_order-i_order_nnt==1 && i_seasonal && i_periodic
            prev_n=n;
        else
            prev_n=0;
        end    
        tobj_seasonal=tab{findtab('SeasonalEffects')};
        data_seasonal=get(findobj(tobj_seasonal,'Tag','uitable_seasonal'),'Data');
        
        line{length(line)+1}=sprintf('  Mode #%i (Seasonal Effects)',i_mode_Seasonal);
        set(edm, 'String',line);
        for n=1:size(data_seasonal,1)
            i_t_pat{n}=str2num(cell2mat(data_seasonal(n,1)));
            i_t_prd{n}=str2num(cell2mat(data_seasonal(n,2)));
            i_t_gr{n}=str2num(cell2mat(data_seasonal(n,3)));
            line{length(line)+1}=sprintf('   Seasonal factor #%i with pattern of length "%i" , growth rate of "%f" and noise "%f"',n,length(i_t_prd{n}),i_t_gr{n},i_ns );
            set(edm, 'String',line);
            [~,U0{i_mode_Seasonal}(:,n+prev_n)]=SimPeriodPattern(i_I_Seasonal,i_t_pat{n},i_t_prd{n},i_ns,i_t_gr{n});
        end

    end
    
    
    % Non-negativity
    tobj_nonnegative=tab{findtab('NonNegative')};
    nng_factors=get(findobj(tobj_nonnegative,'Tag','radiobutton_factors'),'Value');
    nng_whole=get(findobj(tobj_nonnegative,'Tag','radiobutton_whole'),'Value');
    nng_factors_th=str2num(get(findobj(tobj_nonnegative,'Tag','edit_factors'),'String'));
    nng_whole_th=str2num(get(findobj(tobj_nonnegative,'Tag','edit_whole'),'String'));
    
    if nng_factors
        line{length(line)+1}=sprintf(' Applying non-negativity constraint (U{i}<%f=0) on factors ...',nng_factors_th);
        set(edm, 'String',line);
        for n=1:i_order
            U0_nng=U0;
            U0_nng{n}( U0{n} < nng_factors_th ) = 0;
        end 
    end    
           
    
    if i_sparse==1
        line{length(line)+1}=sprintf(' Applying Sparsity on factors...');
        set(edm, 'String',line);
        for i=1:i_order
            line{length(line)+1}=sprintf('  Mode #i: sparsity ratio: %f',i,i_sparse_param);
            set(edm, 'String',line);
            for j=1:size(U0{i},2)
                %nnz1=nnz(U0{i}(:,j));
                [U0{i}(:,j),OS{i}(:,j),OMS{i}(:,j)]=remove_random_elements(U0{i}(:,j),i_sparse_param,0);
                %nnz2=nnz(U0{i}(:,j));
                %line{length(line)+1}=sprintf('  Factor #%i: %f % Sparsed',j,nnz2/nnz1);
                %set(edm, 'String',line);
            end    
            U0{i}=sptensor(U0{i});
        end    
        
    end
    
    if strcmp(i_struct,'Tucker')
        G=reshape(i_lambda,i_R)
        X0=ttensor(tensor(G),U0);
    else
        X0=ktensor(i_lambda',U0);
    end    
    
    if i_normal
            if (isa(X0,'ktensor') || isa(X0,'ttensor')) && i_sparse==0
                X0 = arrange(X0);
                line{length(line)+1}=sprintf(' Normalizing Final Tensor ...');
                set(edm, 'String',line);
            else
                line{length(line)+1}=sprintf(' Warning! Normalizing operation was aborted.');
                set(edm, 'String',line);
            end    
    end    

    if i_fixsign
            if (isa(X0,'ktensor') || isa(X0,'ttensor')) && i_sparse==0
                line{length(line)+1}=sprintf(' Fixing Signs of Final Tensor ...');
                set(edm, 'String',line);
                X0 = fixsigns(X0);
            else
                line{length(line)+1}=sprintf(' Warning! Fixing Sign operation was aborted.');
                set(edm, 'String',line);
            end    
    end

   if i_sparse==0 || i_sparse==2
        X0=full(X0);
   end 
   
   if i_streaming
        if i_temporal==0 || i_nonnegative==0
            line{length(line)+1}=sprintf(' Applying Streaming strucutre with variation parameter of "%f"',i_streaming);
            set(edm, 'String',line);
            X0=SimStream(U0,i_I,i_R,i_streaming);
        else
            line{length(line)+1}=sprintf(' Warning! Streaming strucutre was ignored due to existance of temporal modes');
            set(edm, 'String',line);
        end    
    end
    
    if i_defiend_changes
    if ~isa(U0{1},'sptensor')

      
        tobj_changes=tab{findtab('ChangePoints')};
        data_changes_p=get(findobj(tobj_changes,'Tag','uitable_changepoints_p'),'Data');
        data_changes_s=get(findobj(tobj_changes,'Tag','uitable_changepoints_s'),'Data');
        
        changes_onfactor=get(findobj(tobj_changes,'Tag','changes_onfactor'),'Value');
        changes_ontensor=get(findobj(tobj_changes,'Tag','changes_ontensor'),'Value');
        
        if changes_onfactor
            cmde='factors';
        else
            cmde='final tensor';
        end
                                
        line{length(line)+1}=sprintf(' Adding %i Change-points on %s ...',i_defiend_changes,cmde);
        set(edm, 'String',line);

        for n=1:size(data_changes_p,1)
            i_h_cp_from{n}=cell2mat(data_changes_p(n,1));
            i_h_cp_to{n}=cell2mat(data_changes_p(n,2));
            i_h_cp_val{n}=cell2mat(data_changes_p(n,4));
            i_h_cp_typ{n}=cell2mat(data_changes_p(n,3));

            mu(n)= mean(U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n}));
            s(n)= std(U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n}));
        
            if strfind(i_h_cp_val{n},'mu')
                i_h_cp_val2{n} = strrep(i_h_cp_val{n}, 'mu', '');
                i_h_cp_val2{n}=str2num(i_h_cp_val2{n})*mu(n);
            elseif strfind(i_h_cp_val{n},'sd')
                i_h_cp_val2{n} = strrep(i_h_cp_val{n}, 'sd', '');
                i_h_cp_val2{n}=str2num(i_h_cp_val2{n})*s(n);
            else
                i_h_cp_val2{n}=str2num(i_h_cp_val{n});
            end
            
            for k=1:i_R_Periodic
               if strcmp(i_h_cp_typ{n},'Multiply')
                    if changes_onfactor
                        U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n},k)=U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n},k)*i_h_cp_val2{n};
                    else
                        idx = repmat({':'}, 1, length(size(X0)));
                        idx{i_mode_Periodic}=i_h_cp_from{n}:i_h_cp_to{n};
                        X0(idx{:})=X0(idx{:})*i_h_cp_val2{n};
                    end    
                    if k==1
                        line{length(line)+1}=sprintf('   Mode #%i (Periodic waves): Multiply by %f (%s) the elements in the period between %i to %i',i_mode_Periodic,i_h_cp_val2{n},i_h_cp_val{n},i_h_cp_from{n},i_h_cp_to{n});
                        set(edm, 'String',line);
                    end    
               else
                    if changes_onfactor
                        U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n},k)=U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n},k)+i_h_cp_val2{n};
                    else
                        idx = repmat({':'}, 1, length(size(X0)));
                        idx{i_mode_Periodic}=i_h_cp_from{n}:i_h_cp_to{n};
                        X0(idx{:})=X0(idx{:})+i_h_cp_val2{n};
                    end    
                    if k==1
                        line{length(line)+1}=sprintf('   Mode #%i (Periodic waves): Add %f (%s) to the elements in the period between %i to %i',i_mode_Periodic,i_h_cp_val2{n},i_h_cp_val{n},i_h_cp_from{n},i_h_cp_to{n});
                        set(edm, 'String',line);
                    end
               end
            end
            
        end

       
        for n=1:size(data_changes_s,1)
            i_h_cp_from{n}=cell2mat(data_changes_s(n,1));
            i_h_cp_to{n}=cell2mat(data_changes_s(n,2));
            i_h_cp_val{n}=cell2mat(data_changes_s(n,4));
            i_h_cp_typ{n}=cell2mat(data_changes_s(n,3));

            mu(n)= mean(U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n}));
            s(n)= std(U0{i_mode_Periodic}(i_h_cp_from{n}:i_h_cp_to{n}));
        
            if strfind(i_h_cp_val{n},'mu')
                i_h_cp_val2{n} = strrep(i_h_cp_val{n}, 'mu', '');
                i_h_cp_val2{n}=str2num(i_h_cp_val2{n})*mu(n);
            elseif strfind(i_h_cp_val{n},'sd')
                i_h_cp_val2{n} = strrep(i_h_cp_val{n}, 'sd', '');
                i_h_cp_val2{n}=str2num(i_h_cp_val2{n})*s(n);
            else
                i_h_cp_val2{n}=str2num(i_h_cp_val{n});
            end
            
        
            for k=1:i_R_Seasonal
               if strcmp(i_h_cp_typ{n},'Multiply')
                    if changes_onfactor
                        U0{i_mode_Seasonal}(i_h_cp_from{n}:i_h_cp_to{n},k)=U0{i_mode_Seasonal}(i_h_cp_from{n}:i_h_cp_to{n},k)*i_h_cp_val2{n};
                    else
                        idx = repmat({':'}, 1, length(size(X0)));
                        idx{i_mode_Periodic}=i_h_cp_from{n}:i_h_cp_to{n};
                        X0(idx{:})=X0(idx{:})*i_h_cp_val2{n};
                    end    
                    if k==1
                        line{length(line)+1}=sprintf('   Mode #%i (Seasonal effects): Multiply by %f (%s) the elements in the period between %i to %i',i_mode_Seasonal,i_h_cp_val2{n},i_h_cp_val{n},i_h_cp_from{n},i_h_cp_to{n});
                        set(edm, 'String',line);
                    end    
               else
                    if changes_onfactor
                        U0{i_mode_Seasonal}(i_h_cp_from{n}:i_h_cp_to{n},k)=U0{i_mode_Seasonal}(i_h_cp_from{n}:i_h_cp_to{n},k)+i_h_cp_val2{n};
                    else
                        idx = repmat({':'}, 1, length(size(X0)));
                        idx{i_mode_Periodic}=i_h_cp_from{n}:i_h_cp_to{n};
                        X0(idx{:})=X0(idx{:})+i_h_cp_val2{n};
                    end    
                    if k==1
                        line{length(line)+1}=sprintf('   Mode #%i (Seasonal effects): Add %f (%s) to the elements in the period between %i to %i',i_mode_Seasonal,i_h_cp_val2{n},i_h_cp_val{n},i_h_cp_from{n},i_h_cp_to{n});
                        set(edm, 'String',line);
                    end    
               end
            end
            
        end
        
    else
        line{length(line)+1}=sprintf(' Warning! Change-points were aborted ...');
        set(edm, 'String',line);
    end
    
    end
    
    if i_defiend_anomalies
    if ~isa(U0{1},'sptensor')
        
       tobj_anomalies=tab{findtab('Anomalies')};
       data_anomalies=get(findobj(tobj_anomalies,'Tag','uitable_anomalies'),'Data');
       anomalies_onfactors=get(findobj(tobj_anomalies,'Tag','anomalies_onfactors'),'Value');
       anomalies_ontensor=get(findobj(tobj_anomalies,'Tag','anomalies_ontensor'),'Value');

           
       if anomalies_onfactors
            cmde='factors';
       else
            cmde='final tensor';
       end
        
       line{length(line)+1}=sprintf(' Adding "%i" anomalies on %s ...',i_defiend_changes,cmde);
       set(edm, 'String',line);
        
       for n=1:size(data_anomalies,1)
           anomaly_period_start=str2num(cell2mat(data_anomalies(n,1)));
           anomaly_period_end=str2num(cell2mat(data_anomalies(n,2)));
           tensor_random_generator=strtrim(char(cellstr(data_anomalies(n,3))));
           AnomalyTensorR=str2num(cell2mat(data_anomalies(n,4)));
           if length(AnomalyTensorR)==1
                AnomalyTensorRs=repmat(AnomalyTensorR,length(anomaly_period_start),1);
           end  
           anomaly_lambda_factor=cell2mat(data_anomalies(n,5));
           anomaly_lambda=i_lambda*anomaly_lambda_factor;
           AnomalyTensorSize=anomaly_period_end-anomaly_period_start+1;
           line{length(line)+1}=sprintf('  Anomaly #%i: Random generator: %s, size: [%s], rank: %i, lambda/core: %s',n,tensor_random_generator,strrep(num2str(AnomalyTensorSize),'  ',' '),AnomalyTensorR,strrep(num2str(anomaly_lambda),'  ',' '));
           set(edm, 'String',line);
           
           for i=1:length(AnomalyTensorSize)
               set(edm, 'String',line);
               AnomalyTensorU{i}=SimFacGen(tensor_random_generator,AnomalyTensorSize(i),AnomalyTensorRs(i),0,0);
               if anomalies_onfactors
                  U0{i}(anomaly_period_start(i):anomaly_period_end(i),:)=AnomalyTensorU{i};
               end
           end

           if anomalies_ontensor
                X_anomaly=ktensor(anomaly_lambda',AnomalyTensorU);
                for i=1:size(anomaly_period_start,2)
                    anomaly_period(i)={anomaly_period_start(i):anomaly_period_end(i)};
                end
                X0(anomaly_period{:})=X_anomaly;
           end           

       end
       
    else
        line{length(line)+1}=sprintf(' Warning! Anomalies were aborted ...');
        set(edm, 'String',line);
    end
    end

    X=X0;
    
    % Noise
    tobj=tab{findtab('Noise')};
    radiobutton_whole=get(findobj(tobj,'Tag','radiobutton_whole'),'Value');
    edit_noise_snr=str2num(get(findobj(tobj,'Tag','edit_noise_snr'),'String'));
    radiobutton_sparse=get(findobj(tobj,'Tag','radiobutton_sparse'),'Value');
    edit_noise_level=str2num(get(findobj(tobj,'Tag','edit_noise_level'),'String'));
    
    if radiobutton_whole==1
         line{length(line)+1}=sprintf(' Adding white noise (%f dB) on final tensor ...',edit_noise_snr);
         set(edm, 'String',line);
         X=add_awgn_noise(X0,edit_noise_snr); % by Mathuranathan Viswanathan 
    end    
    if radiobutton_sparse==1
         line{length(line)+1}=sprintf(' Adding sparse white noise (ratio of %f) on final tensor ...',edit_noise_level);
         set(edm, 'String',line);
         X=add_sparse_noise(X0,edit_noise_level,1);
    end 

    % Non-negative
    if nng_whole
        line{length(line)+1}=sprintf(' Applying non-negativity constraint (X<%f=0) on final tensor ...',nng_whole_th);
        set(edm, 'String',line);
        X_nng=reshape(double(X),prod(size(X)),1);
        X_nng(X_nng<nng_whole_th)=0;
        X=tensor(reshape(X_nng,size(X)));
    end
    
    % Sparsity
    tobj_sparsity=tab{findtab('Sparsity')};
    radiobutton_sparsetensor=get(findobj(tobj_sparsity,'Tag','radiobutton_sparsetensor'),'Value');
    i_m_rate=str2num(get(findobj(tobj_sparsity,'Tag','edit_sparsetensor'),'String'));
    checkbox_sparsetensor=get(findobj(tobj_sparsity,'Tag','checkbox_sparsetensor'),'Value');
    radiobutton_counttensor=get(findobj(tobj_sparsity,'Tag','radiobutton_counttensor'),'Value');
    sparsetensor_R=str2num(get(findobj(tobj_sparsity,'Tag','sparsetensor_R'),'String'));
    sparsetensor_R2=str2num(get(findobj(tobj_sparsity,'Tag','sparsetensor_R2'),'String'));
    sparsetensor_ceiling=str2num(get(findobj(tobj_sparsity,'Tag','sparsetensor_ceiling'),'String'));
   
    if radiobutton_sparsetensor
        line{length(line)+1}=sprintf(' Applying sparsity on final tensor (ratio: %f) ...',i_m_rate);
        set(edm, 'String',line);
        [X,O,Omega]=remove_random_elements(X,i_m_rate,checkbox_sparsetensor);
    end    
    
    if radiobutton_counttensor
        line(2:length(line))=[];
        line{length(line)+1}=sprintf('Warning! All confirguations were ignored, because you chose to generate sparse count tensor (Sparsity tab)',num2str(i_I),sparsetensor_R,sparsetensor_R2);
        line{length(line)+1}=sprintf('Generating sparse count tensor (Size: [%s], Rank: %i, Noise components: %i)',num2str(i_I),sparsetensor_R,sparsetensor_R2);
        set(edm, 'String',line);
        lambda=SimLamGen('CP',sparsetensor_R,'gamma',100,1);
        lambda(sparsetensor_R+1:sparsetensor_R2+sparsetensor_R+1-1)=rand(sparsetensor_R2,1);
        for n=1:length(i_I);
            [~, U0{n}]=SimFacGen('gamma',i_I(n),sparsetensor_R+sparsetensor_R2,0,1);
        end
        X0=ktensor(lambda',U0);
        X=sptensor(poissrnd(double(X0))); 
        U0=[];
        X0=[];
    end    

    line{length(line)+1}=sprintf('The tensor was generated sucessfully! \r\n');
    set(edm, 'String',line);
    set(savebutton,'Enable','On');
        
    
end    
    
function generate_edit(hObject, eventdata)
    eventdata.PreviousData
     input = get(hObject,'String');
     set(eventdata.Source,'String',EventData.PreviousData);
end
    
function  sparsity_update(PushButton, EventData)
    tobj=tab{findtab('Sparsity')};
    radiobutton_factors=get(findobj(tobj,'Tag','radiobutton_factors'),'Value');
    radiobutton_sparsetensor=get(findobj(tobj,'Tag','radiobutton_sparsetensor'),'Value');
    radiobutton_counttensor=get(findobj(tobj,'Tag','radiobutton_counttensor'),'Value');
    edit_factors=findobj(tobj,'Tag','edit_factors');
    edit_sparsetensor=findobj(tobj,'Tag','edit_sparsetensor');
    checkbox_sparsetensor=findobj(tobj,'Tag','checkbox_sparsetensor');
    sparsetensor_R=findobj(tobj,'Tag','sparsetensor_R');
    sparsetensor_R2=findobj(tobj,'Tag','sparsetensor_R2');
    sparsetensor_ceiling=findobj(tobj,'Tag','sparsetensor_ceiling');
    set(edit_factors,'Enable',offon(radiobutton_factors));
    set(edit_sparsetensor,'Enable',offon(radiobutton_sparsetensor));
    set(checkbox_sparsetensor,'Enable',offon(radiobutton_sparsetensor));
    set(sparsetensor_R,'Enable',offon(radiobutton_counttensor));
    set(sparsetensor_R2,'Enable',offon(radiobutton_counttensor));
    set(sparsetensor_ceiling,'Enable',offon(radiobutton_counttensor));
    if radiobutton_factors==1
        i_sparse=1;
        i_sparse_param=str2num(get(edit_factors,'String'));
    elseif  radiobutton_sparsetensor==1
        i_sparse=2;
    elseif  radiobutton_counttensor==1
        i_sparse=3;
        i_sparse_param=edit_factors;
    else
        i_sparse=0;
    end
    status_update();
end

function nonnegative_update(PushButton, EventData)
	tobj=tab{findtab('NonNegative')};
    radiobutton_none=get(findobj(tobj,'Tag','radiobutton_none'),'Value');
    radiobutton_factors=get(findobj(tobj,'Tag','radiobutton_factors'),'Value');
    radiobutton_whole=get(findobj(tobj,'Tag','radiobutton_whole'),'Value');
    edit_factors=findobj(tobj,'Tag','edit_factors');
    edit_whole=findobj(tobj,'Tag','edit_whole');
    set(edit_factors,'Enable',offon(radiobutton_factors));
    set(edit_whole,'Enable',offon(radiobutton_whole));
    if radiobutton_factors==1 || radiobutton_whole==1 
        i_nonnegative=1;
    else
        i_nonnegative=0;
    end
    status_update();
end 

function noise_update(PushButton, EventData)
    tobj=tab{findtab('Noise')};
    
    radiobutton_nonoise=get(findobj(tobj,'Tag','radiobutton_nonoise'),'Value');
    radiobutton_factors=get(findobj(tobj,'Tag','radiobutton_factors'),'Value');
    radiobutton_whole=get(findobj(tobj,'Tag','radiobutton_whole'),'Value');
    radiobutton_sparse=get(findobj(tobj,'Tag','radiobutton_sparse'),'Value');
    edit_noise_factors=findobj(tobj,'Tag','edit_noise_factors');
    edit_noise_snr=findobj(tobj,'Tag','edit_noise_snr');
    edit_noise_level=findobj(tobj,'Tag','edit_noise_level');

    set(edit_noise_factors,'Enable',offon(radiobutton_factors));
    set(edit_noise_snr,'Enable',offon(radiobutton_whole));
    set(edit_noise_level,'Enable',offon(radiobutton_sparse));
    
    if radiobutton_factors==1 || radiobutton_whole==1 || radiobutton_sparse==1
        i_noise=1;
    else
        i_noise=0;
    end
    
    if radiobutton_factors==1
        i_ns=str2num(get(edit_noise_factors,'String'));
    else
        i_ns=0;
    end    
    
    status_update();
    
end

function anomalies_empty(PushButton, EventData)
    tobj=tab{findtab('Anomalies')};
    uitable_anomalies=findobj(tobj,'Tag','uitable_anomalies');
    set(uitable_anomalies,'Data',{blanks(0),blanks(0),blanks(0),blanks(0),blanks(0)});
    set(uitable_anomalies,'ColumnEditable',logical([0 0 0 0 0]) );
    i_defiend_anomalies=0;
    status_update();
end
function anomalies_edit(PushButton, EventData)
    returnback=0;
    col=EventData.Indices(2);
    rw=EventData.Indices(1);
    dat=get(EventData.Source,'Data');

    if col==1 || col==2
        prv=str2num(EventData.PreviousData);
        nwe=str2num(EventData.NewData);
        col1=str2num(dat{rw,1});
        col2=str2num(dat{rw,2});
        if isnan(nwe)
             returnback=1;
        end
        if length(nwe)~= length(i_I)
             returnback=1;
        end 
        if ~isempty(find(nwe<1))
            returnback=1;
        end 
        if returnback==0
            if ~isempty(find(i_I-nwe <0))
                returnback=1;
            end    
            if ~isempty(find(col2-col1<0))
                returnback=1;
            end 
        end    
    end
        
    if col==4
        prv=str2num(EventData.PreviousData);
        nwe=str2num(EventData.NewData);
        if isempty(nwe)
            returnback=1;
        end    
        if length(nwe)>1 && length(nwe)~=length(i_I)
             returnback=1;
        end    
        if ~isempty(find(nwe<1))
              returnback=1;
        end    
       
    end
    if col==5
       prv=EventData.PreviousData;
       nwe=EventData.NewData;
       if isnan(nwe)
            returnback=1;
       end
    end
    
     if returnback==1       
        dat{rw,col}=num2str(prv);
        set(EventData.Source,'Data',dat);
        errordlg('Invalid edit! the value was changed to its previous state');
     end
end

function anomalies_slide(PushButton, EventData)
    tobj=tab{findtab('Anomalies')};
    anomalies_length=findobj(tobj,'Tag','anomalies_length');
    anomalies_slide_prev=findobj(tobj,'Tag','anomalies_slide_prev');
    anomalies_length=get(anomalies_length,'Value');
    small_tensor_size=ceil(i_I*anomalies_length);
    anomalies_slide_prev.String=mat2str(small_tensor_size);
end
function anomalies_add(PushButton, EventData)
    tobj=tab{findtab('Anomalies')};
    tbl=findobj(tobj,'Tag','uitable_anomalies');
    anomalies_length=get(findobj(tobj,'Tag','anomalies_length'),'Value');
    lambdafactor=str2num(get(findobj(tobj,'Tag','lambdafactor'),'String'));
    anomalies_factorgenerator=get(findobj(tobj,'Tag','anomalies_factorgenerator'),'Value');
    obj3=findobj(tobj,'Tag','anomalies_factorgenerator');
    fac_methods=cellstr(get(obj3,'String'))';
    anomalies_tensorrank=str2num(get(findobj(tobj,'Tag','anomalies_tensorrank'),'String'));
    anomalies_length=ceil(i_I*anomalies_length);
    dat=get(tbl,'Data');
    if all(cellfun(@isempty, dat(:)))    
       dat=[];
       set(tbl,'Data',dat);
       rng_max=i_I;
    else
       rng_max=str2num(dat{size(dat,1),1})-1;
    end   
    if ~isempty(find(rng_max-anomalies_length+1<1))
       errordlg('It is not allowed to add more anomalies!');
       return;
    end
    
    row{1,1}=strrep(strrep(mat2str(rng_max-anomalies_length+1),'[',''),']','');
    row{1,2}=strrep(strrep(mat2str(rng_max),'[',''),']','');
    row{1,3}=fac_methods{anomalies_factorgenerator};
    row{1,4}=mat2str(anomalies_tensorrank);
    row{1,5}=lambdafactor;
    dat=[dat; row];
    set(tbl,'ColumnEditable',logical([1 1 1 1 1]) );
    set(tbl, 'ColumnFormat', {'numeric' 'numeric' fac_methods 'numeric' 'char'});
    set(tbl,'Data',dat);
     i_defiend_anomalies=i_defiend_anomalies+1;
    status_update();
end

function changes_edit(PushButton, EventData)
    returnback=0;
    col=EventData.Indices(2);
    if col==4
        rw=EventData.Indices(1);
        prv=EventData.PreviousData;
        nwe=EventData.NewData;
        if strfind(nwe,'mu')
            newn2=sprintf('%0.2fmu',str2num(strrep(nwe, 'mu', '')));
        elseif strfind(nwe,'sd')  
            newn2=sprintf('%0.2fsd',str2num(strrep(nwe, 'sd', '')));
        elseif strfind(nwe,'i')  
            newn2=sprintf('%0.2f',str2num(strrep(nwe, 'i', '')));
        else
            newn2=sprintf('%0.2f',str2num(nwe));
            if isempty(str2num(nwe))
               newn2=EventData.PreviousData;
               errordlg('Invalid edit! the value was changed to its previous state');
            end
        end  
        newn2=strrep(newn2, ' ', '');
        dat=get(EventData.Source,'Data');
        dat(rw,col)={newn2};
        set(EventData.Source,'Data',dat);
    elseif col==1 || col==2
        rw=EventData.Indices(1);
        prv=EventData.PreviousData;
        nwe=EventData.NewData;
        if strfind(PushButton.Tag,'_p')
            mx=i_I_Periodic;
        else
            mx=i_I_Seasonal;
        end    
        
        dat=get(EventData.Source,'Data');
        if nwe<1 || nwe>mx
            dat(rw,col)={prv};
            errordlg(sprintf('The value should be between 1 and %i! it was changed to its previous state',mx));
            set(EventData.Source,'Data',dat);
        end
        
        if cell2mat(dat(rw,1)) >= cell2mat(dat(rw,2))
          errordlg(sprintf('The start point should be lower than the end point! The value was changed to its previous state',mx));
          dat(rw,col)={prv};
          set(EventData.Source,'Data',dat);
        end    
        
    end

    
end    

function changes_reset(PushButton, EventData)
    tobj=tab{findtab('ChangePoints')};
    tbl_p=findobj(tobj,'Tag','uitable_changepoints_p');
    tbl_s=findobj(tobj,'Tag','uitable_changepoints_s');
    set(tbl_s,'Data',{blanks(0),blanks(0),blanks(0),blanks(0)});
    set(tbl_p,'Data',{blanks(0),blanks(0),blanks(0),blanks(0)});
    i_defiend_changes=0;
    status_update();
end

function changes_add(PushButton, EventData)
    tobj=tab{findtab('ChangePoints')};
    cp_number_of_changes=str2num(get(findobj(tobj,'Tag','cp_number_of_changes'),'String'));
    cp_seasonal=get(findobj(tobj,'Tag','cp_seasonal'),'Value');
    cp_periodic=get(findobj(tobj,'Tag','cp_periodic'),'Value');
    cp_length=get(findobj(tobj,'Tag','cp_length'),'Value');
    cp_multiply=get(findobj(tobj,'Tag','cp_multiply'),'Value');
    cp_add=get(findobj(tobj,'Tag','cp_add'),'Value');
    cp_mix=get(findobj(tobj,'Tag','cp_mix'),'Value');
    cp_typ(1)=get(findobj(tobj,'Tag','cp_constant'),'Value');
    cp_typ(2)=get(findobj(tobj,'Tag','cp_sd_factor'),'Value');
    cp_typ(3)=get(findobj(tobj,'Tag','cp_mu_factor'),'Value');
    cp_constant_from=str2num(get(findobj(tobj,'Tag','cp_constant_from'),'String'));
    cp_constant_to=str2num(get(findobj(tobj,'Tag','cp_constant_to'),'String'));
    cp_sd_from=str2num(get(findobj(tobj,'Tag','cp_sd_from'),'String'));
    cp_sd_to=str2num(get(findobj(tobj,'Tag','cp_sd_to'),'String'));
    cp_mu_from=str2num(get(findobj(tobj,'Tag','cp_mu_from'),'String'));
    cp_mu_to=str2num(get(findobj(tobj,'Tag','cp_mu_to'),'String'));
    cp_constant_int=get(findobj(tobj,'Tag','cp_constant_int'),'Value');
    cp_sd_int=get(findobj(tobj,'Tag','cp_sd_int'),'Value');
    cp_mu_int=get(findobj(tobj,'Tag','cp_mu_int'),'Value');

    if cp_multiply==1
        for i=1:cp_number_of_changes
            cp_operator{i}='Multiply';
        end
    elseif cp_add==1;
        for i=1:cp_number_of_changes
           cp_operator{i}='Add';
        end   
    else
        for i=1:cp_number_of_changes
            if randi([1 2],1,1)==1
                cp_operator{i}='Add';
            else
                cp_operator{i}='Multiply';
            end
        end
    end    
    
    if cp_constant_int==0
        fparams{1}=(cp_constant_to-cp_constant_from).*rand(1,cp_number_of_changes)+cp_constant_from;
    elseif  cp_constant_int==1
        fparams{1}=randi([cp_constant_from cp_constant_to],1,cp_number_of_changes);
    end    
    if cp_sd_int==0
        fparams{2}=(cp_sd_to-cp_sd_from).*rand(1,cp_number_of_changes)+cp_sd_from;
    elseif  cp_sd_int==1
        fparams{2}=randi([cp_sd_from cp_sd_to],1,cp_number_of_changes);
    end    
    
    if cp_mu_int==0
        fparams{3}=(cp_mu_to-cp_mu_from).*rand(1,cp_number_of_changes)+cp_mu_from;
    elseif  cp_mu_int==1
        fparams{3}=randi([cp_mu_from cp_mu_to],1,cp_number_of_changes);
    end  
    
    tidx=find(cp_typ==1);
    rts=randi([min(tidx) max(tidx)],1,cp_number_of_changes);
    for i=1:cp_number_of_changes
        pms{i}=fparams{rts(i)}(i);
        if rts(i)==1
            pms{i}=sprintf('%0.2f',pms{i});
        elseif rts(i)==2
            pms{i}=sprintf('%0.2fsd',pms{i});
        elseif rts(i)==3
            pms{i}=sprintf('%0.2fmu',pms{i});
           end    
    end
    tbl_p=findobj(tobj,'Tag','uitable_changepoints_p');
    tbl_s=findobj(tobj,'Tag','uitable_changepoints_s');
    changepoints_p_title=findobj(tobj,'Tag','changepoints_p_title');
    changepoints_s_title=findobj(tobj,'Tag','changepoints_s_title');
  
    cp_length;
    cp_length_s=ceil(i_I_Seasonal*cp_length);
    cp_length_p=ceil(i_I_Periodic*cp_length);  

    if cp_seasonal==1
        startpoints_s=randperiods(i_I_Seasonal,cp_length_s,cp_number_of_changes);
        endpoints_s=startpoints_s+cp_length_s;
    end
    if cp_periodic==1
        startpoints_p=randperiods(i_I_Periodic,cp_length_p,cp_number_of_changes);
        endpoints_p=startpoints_p+cp_length_p;
    end    
    
    if cp_seasonal==1
        for i=1:cp_number_of_changes
            if endpoints_s(i)>i_I_Seasonal
                endpoints_s(i)=0;
                startpoints_s(i)=0;
            end    
            if startpoints_s(i)==0
                startpoints_s(i)=0;
                endpoints_s(i)=0;
            end    
        end    
        startpoints_s(find(startpoints_s==0))=[];
        endpoints_s(find(endpoints_s==0))=[];
    end        
    
    if cp_periodic==1
        for i=1:cp_number_of_changes
            if endpoints_p(i)>i_I_Periodic
                endpoints_p(i)=0;
                startpoints_p(i)=0;
            end 
            if startpoints_p(i)==0
                startpoints_p(i)=0;
                endpoints_p(i)=0;
            end 
        end 
        startpoints_p(find(startpoints_p==0))=[];
        endpoints_p(find(endpoints_p==0))=[];
    end


    dat_p={blanks(0),blanks(0),blanks(0),blanks(0)};
    if cp_periodic==1
        set(tbl_p,'Visible','on');
        set(changepoints_p_title,'Visible','on');
        for i=1:length(startpoints_p)
            dat_p(i,1)={startpoints_p(i)};
            dat_p(i,2)={endpoints_p(i)};
            dat_p(i,3)={cp_operator{i}};
            dat_p(i,4)={pms{i}};
        end
    end
    
    dat_s={blanks(0),blanks(0),blanks(0),blanks(0)};
    if cp_seasonal==1
        set(tbl_s,'Visible','on');
        set(changepoints_s_title,'Visible','on');
        for i=1:length(startpoints_s)
            dat_s(i,1)={startpoints_s(i)};
            dat_s(i,2)={endpoints_s(i)};
            dat_s(i,3)={cp_operator{i}};
            dat_s(i,4)={pms{i}};
        end
    end

    set(tbl_s,'Data',dat_s);
    set(tbl_p,'Data',dat_p);

    
    for i=1:cp_number_of_changes
        RowsN_p(i,1)={['ChangePoint ',sprintf('%i',i)]};
    end    
    set(tbl_p,'RowName',RowsN_p);
    
    for i=1:cp_number_of_changes
        RowsN_s(i,1)={['ChangePoint ',sprintf('%i',i)]};
    end  
    
    set(tbl_s,'RowName',RowsN_s);
    
    set(tbl_s,'ColumnEditable',logical([1 1 1 1]) );
    set(tbl_p,'ColumnEditable',logical([1 1 1 1]) );
    if cp_seasonal==0
       set(tbl_s,'Visible','off');
       set(changepoints_s_title,'Visible','off');
    end
    if cp_periodic==0
       set(tbl_p,'Visible','off');
       set(changepoints_p_title,'Visible','off');
    end
    i_defiend_changes=cp_number_of_changes;

    status_update();
    
end    

function sp=randperiods(i_I_Periodic,cp_length_p,cp_number_of_changes)
    sp(1)=ceil(i_I_Periodic/2);
    sp(2)=1;
    sp(3)=i_I_Periodic-cp_length_p;
    sp(4)=ceil((sp(1)-sp(2))/2);
    sp(5)=ceil(sp(3)-ceil((sp(3)-sp(1))/2));
    sp(6)=ceil((sp(1)-sp(4))/2);
    sp(7)=ceil(sp(3)-ceil((sp(3)-sp(5))/2));
    for i=8:2:cp_number_of_changes
       sp(i)=ceil((sp(1)-sp(i-2))/2);
       sp(i+1)=ceil(sp(3)-ceil((sp(3)-sp(i-1))/2));
    end
    sp=sp(1:cp_number_of_changes);
end    
    
function out=offon(val)
    if isempty(val)
        val=0;
    end    
    if val==0
        out='off';
    else
       out='on';         
    end
end        

function TabNo=findtab(str)
    TabNo=1;
    for i = 1 : length(tabs)
        if strcmp(str,tabs{i})
            TabNo=i;
            break;
            %return;
        end    
    end
end 

function status_update()
 
    set(findobj(sb,'Tag','Strcutrue'),'String',['Strcutrue: ',num2str(i_struct)]);
    set(findobj(sb,'Tag','Size'),'String',['Size: ',num2str(i_I)]);
    set(findobj(sb,'Tag','Rank'),'String',['Rank: ',num2str(i_R)]);
    set(findobj(sb,'Tag','Tensor Order'),'String',['Tensor Order: ',num2str(i_order)]);
    set(findobj(sb,'Tag','Temporal'),'Enable',offon(i_temporal));
    set(findobj(sb,'Tag','Periodic'),'Enable',offon(i_defiend_periodic));
    set(findobj(sb,'Tag','Seasonal'),'Enable',offon(i_defiend_seasonal));
    set(findobj(sb,'Tag','Changes'),'Enable',offon(i_defiend_changes));
    set(findobj(sb,'Tag','Anomalies'),'Enable',offon(i_defiend_anomalies));
    set(findobj(sb,'Tag','Fixed Angle'),'Enable',offon(i_s_angle));
    set(findobj(sb,'Tag','Fixed correlation'),'Enable',offon(i_s_corr));
    set(findobj(sb,'Tag','Noise'),'Enable',offon(i_noise));
    set(findobj(sb,'Tag','Non-negative'),'Enable',offon(i_nonnegative));
    set(findobj(sb,'Tag','Sparse'),'Enable',offon(i_sparse));
    set(findobj(sb,'Tag','Normalize'),'Enable',offon(i_normal));
    set(findobj(sb,'Tag','Sign Fixing'),'Enable',offon(i_fixsign));
    if i_s_angle
        set(findobj(sb,'Tag','Fixed Angle'),'String',['Fixed Angle (',num2str(i_s_angle),'°)']);
    else
        set(findobj(sb,'Tag','Fixed Angle'),'String','Fixed Angle');
    end
    if i_s_corr
        set(findobj(sb,'Tag','Fixed correlation'),'String',['Fixed correlation (',num2str(i_s_corr),')']);
    else
        set(findobj(sb,'Tag','Fixed correlation'),'String','Fixed correlation');        
    end    
    if i_defiend_periodic
        set(findobj(sb,'Tag','Periodic'),'String',['Periodic (',num2str(i_defiend_periodic),')']);
    else
        set(findobj(sb,'Tag','Periodic'),'String','Periodic');
    end    
	if i_defiend_seasonal
        set(findobj(sb,'Tag','Seasonal'),'String',['Seasonal (',num2str(i_defiend_seasonal),')']);
    else
        set(findobj(sb,'Tag','Seasonal'),'String','Seasonal');
    end 
    if i_defiend_changes
        set(findobj(sb,'Tag','Changes'),'String',['Changes (',num2str(i_defiend_changes),')']);
    else
        set(findobj(sb,'Tag','Changes'),'String','Changes');
    end
    if i_defiend_anomalies
        set(findobj(sb,'Tag','Anomalies'),'String',['Anomalies (',num2str(i_defiend_anomalies),')']);
    else
        set(findobj(sb,'Tag','Anomalies'),'String','Anomalies');
    end
    
    set(findobj(sb,'Tag','Streaming'),'Enable',offon(i_streaming));

end

function factor_update_status(hndl)
    dat=get(hndl,'Data');
    i_I=cell2mat(dat(:,2))';
    i_R=cell2mat(dat(:,3))';
    i_temporal=~isempty(cell2mat(strfind(dat(:,1),'Temporal')));
    i_seasonal=~isempty(cell2mat(strfind(dat(:,1),'Seasonal')));
    i_periodic=~isempty(cell2mat(strfind(dat(:,1),'Periodic')));
    ip=strfind(dat(:,1),'Periodic');
    is=strfind(dat(:,1),'Seasonal');
    i_mode_Periodic = find(not(cellfun('isempty', ip)));
    i_mode_Seasonal = find(not(cellfun('isempty', is)));
    
    i_I_Seasonal=cell2mat(dat(i_mode_Seasonal,2));
    i_R_Seasonal=cell2mat(dat(i_mode_Seasonal,3));
    i_I_Periodic=cell2mat(dat(i_mode_Periodic,2));
    i_R_Periodic=cell2mat(dat(i_mode_Periodic,3));

    set(findobj(tab{findtab('Modes')},'Tag','tensor_struct_cp'),'Enable','on');
    if length(unique(i_R))>1
        set(findobj(tab{findtab('Modes')},'Tag','tensor_struct_cp'),'Enable','off');
    end
    
    if  get(findobj(tab{findtab('Modes')},'Tag','tensor_struct_tucker'),'Value')==0
    
        if length(unique(i_R))==1
            i_struct='CP';
        else
            i_struct='Tucker';
            set(findobj(tab{findtab('Modes')},'Tag','tensor_struct_cp'),'Enable','off');
            set(findobj(tab{findtab('Modes')},'Tag','tensor_struct_auto'),'Value',1);
        end    
    else
        i_struct='Tucker';
    end
    i_order=size(dat,1);  
    
    status_update();
    modes_update_lambda();

end

function modes_update_lambda (PushButton, EventData)
    % 1 Ones      2 Rand      3 Randn     4 Customized
    if i_added==0
        return;
    end    
	tobj=tab{findtab('Modes')};
    modes_lambda_method=get(findobj(tobj,'Tag','modes_lambda_method'),'Value');
    modes_lambda_method_topn=get(findobj(tobj,'Tag','modes_lambda_method_topn'),'Value');
    modes_lambda_method_weight=str2num(get(findobj(tobj,'Tag','modes_lambda_method_weight'),'String'));
    lambdaobj=findobj(tobj,'Tag','modes_lambda');
    modes_lambda=get(lambdaobj,'String');
    set(lambdaobj,'Enable','off');
    switch modes_lambda_method
        case 1
            i_lambda=SimLamGen(i_struct,i_R,'ones',modes_lambda_method_topn,modes_lambda_method_weight);
        case 2
            i_lambda=SimLamGen(i_struct,i_R,'rand',modes_lambda_method_topn,modes_lambda_method_weight);
        case 3
            i_lambda=SimLamGen(i_struct,i_R,'randn',modes_lambda_method_topn,modes_lambda_method_weight);
        case 4
            i_lambda=str2num(modes_lambda);
            set(lambdaobj,'Enable','on');
    end
    
    if modes_lambda_method~=4
        set(lambdaobj,'String',num2str(i_lambda));
    end


end
    
function modes_cptucker (PushButton, EventData)
    tobj=tab{findtab('Modes')};
    tensor_struct_auto=get(findobj(tobj,'Tag','tensor_struct_auto'),'Value');
    tensor_struct_cp=get(findobj(tobj,'Tag','tensor_struct_cp'),'Value');
    tensor_struct_tucker=get(findobj(tobj,'Tag','tensor_struct_tucker'),'Value');
    if tensor_struct_tucker
        i_struct='Tucker';
    end     
    if tensor_struct_cp
        i_struct='CP';
    end    
    status_update();
    modes_update_lambda();
end

function modes_celledit(PushButton, EventData)
    returnback=0;
    col=EventData.Indices(2);
    if col==1
        rw=EventData.Indices(1);
        prv=EventData.PreviousData;
        nwe=EventData.NewData;
        prv_nt=isempty(strfind(char(cellstr(prv)),'Non-temporal'));
        new_nt=isempty(strfind(char(cellstr(nwe)),'Non-temporal'));
        if new_nt==1
            returnback=1;
        end
        
        if returnback~=1
            rw=EventData.Indices(1);
            col=EventData.Indices(2);
            nwe=EventData.NewData;
            isnot_temporal=isempty(strfind(char(cellstr(nwe)),'(Temporal)'));
            if isnot_temporal==1
                i_f_method{rw}=strrep(char(cellstr(nwe)),' (Non-temporal)','');
            end   
        end 
    
    else
        prv=EventData.PreviousData;
        nwe=EventData.NewData;
        if isnan(nwe)
            returnback=1;
        else
            if nwe<1
                returnback=1;
            end    
        end
    end

            
	if returnback==1
        rw=EventData.Indices(1);
        dat=get(EventData.Source,'Data');
        dat(rw,col)={prv};
        set(EventData.Source,'Data',dat);
        errordlg('Invalid edit! the value was changed to its previous state');
    end
    
    
    factor_update_status(EventData.Source);
end        

function modes_add_basic(PushButton, EventData)
 	tobj=tab{findtab('Modes')};
    obj1=findobj(tobj,'Tag','tensorsize');
	obj2=findobj(tobj,'Tag','tensorrank');
	obj3=findobj(tobj,'Tag','factorgenerator');
  	obj4=findobj(tobj,'Tag','temporalfactorgenerator');
    tbl=findobj(tobj,'Tag','uitable_factors');
    set(tbl,'Data',{blanks(0),blanks(0),blanks(0)});
    dat=get(tbl,'Data');
    fac_methods1=cellstr(get(obj3,'String'))';
    for i=1:length(fac_methods1)
        fac_methods{i}=[fac_methods1{i} ' (Non-temporal)'];
    end
    fac_methods{length(fac_methods)+1}='Periodic (Temporal)';
    fac_methods{length(fac_methods)+1}='Seasonal (Temporal)';
    fac_methods{length(fac_methods)+1}='Periodic/Seasonal (Temporal)';
    %save ('facmethods.mat','fac_methods');
    %String = sprintf('Data%d|',1:size(Data))
    tensorsize=str2num(get(obj1, 'String'));
    tensorrank=str2num(get(obj2, 'String'));
    tensorranks=[tensorrank repmat(tensorrank,1,length(tensorsize)-length(tensorrank))];
    factorgenerator=get(obj3, 'Value');
    temporalfactorgenerator=get(obj4, 'Value');
    temporal_modes=[2 1 1 1 0];
    cnt_temp=temporal_modes(temporalfactorgenerator);
    i_order_nnt=length(tensorsize)-cnt_temp;
    for i=1:i_order_nnt;
        dat(i,1)={fac_methods{factorgenerator}};
        i_f_method{i}=fac_methods1{factorgenerator};
        dat(i,2)={tensorsize(i)};
        dat(i,3)={tensorranks(i)};
    end
    for i=length(tensorsize)-cnt_temp+1:length(tensorsize);
        dat(i,2)={tensorsize(i)};
        dat(i,3)={tensorranks(i)};
    end
    switch temporalfactorgenerator
        case 1
             dat(length(tensorsize)-cnt_temp+1,1)={'Periodic (Temporal)'};
             dat(length(tensorsize),1)={'Seasonal (Temporal)'};
        case 2
             dat(length(tensorsize)-cnt_temp+1,1)={'Periodic/Seasonal (Temporal)'};
        case 3
             dat(length(tensorsize)-cnt_temp+1,1)={'Periodic (Temporal)'};
             i_defiend_seasonal=0;
             seasonal_reset('','');
        case 4
             dat(length(tensorsize)-cnt_temp+1,1)={'Seasonal (Temporal)'};
             i_defiend_periodic=0;
             periodic_empty('','');
        otherwise
    end
    set(tbl,'Data',dat);
    for i=1:length(tensorsize)
        RowsN(i,1)={['Mode ',sprintf('%i',i)]};
    end    
    set(tbl,'RowName',RowsN);
    
    set(tbl,'ColumnEditable',logical([1 1 1]) );
    set(tbl, 'ColumnFormat', {fac_methods 'numeric' 'numeric' });
    set (tbl,'ColumnWidth',{ 260 50 50 });

    factor_update_status(tbl);
    i_added=1;
end

function seasonal_scale_check(PushButton, EventData)
	tobj=tab{findtab('SeasonalEffects')};
    st=findobj(tobj,'Tag','seasonal_template');
	sw=findobj(tobj,'Tag','seasonal_weights');
	sd=findobj(tobj,'Tag','seasonal_durations');
    seas_scl=findobj(tobj,'Tag','seasonal_scale');
    
    set(sw,'String','');
    set(sd,'String','');
    set(st,'Value',1);
    if get(seas_scl,'Value')==1
        set(st,'Enable','off');
    else
        set(st,'Enable','on');
    end
    
end

function seasonal_update(PushButton, EventData)
        tobj=tab{findtab('SeasonalEffects')};
        st=findobj(tobj,'Tag','seasonal_template');
        sel=get(st,'Value');
        sw=findobj(tobj,'Tag','seasonal_weights');
        sd=findobj(tobj,'Tag','seasonal_durations');
        ss=findobj(tobj,'Tag','seasonal_scale');
        scle=get(ss,'Value');
        % 2)Hour  3)Day  4)Week  5)Month  6)Season
        % 2 MorningSimple(4)  3 MorningAdvanced (8) 4 Day(24) 9 Year(365) 5 Week(7) 6 Weekend/Weekends (2) 7 Month (12) 8 Season(4)                                    
        switch scle
            case 2 %Hour
                scw{2}=1; scw{3}=1; scw{4}=1; scw{9}=24; scw{5}=24; scw{6}=24; scw{7}=24; scw{8}=24;
            case 3 %Day
                scw{2}=0; scw{3}=0; scw{4}=0; scw{9}=1; scw{5}=1; scw{6}=1; scw{7}=1; scw{8}=1;   
            case 4 %Week
                scw{2}=0; scw{3}=0; scw{4}=0; scw{9}=1; scw{5}=0; scw{6}=0; scw{7}=1/7; scw{8}=1/7;       
            case 5 %Month
                scw{2}=0; scw{3}=0; scw{4}=0; scw{9}=1; scw{5}=0; scw{6}=0; scw{7}=0; scw{8}=1/28;       
            case 6 %Season
                scw{2}=0; scw{3}=0; scw{4}=0; scw{9}=1;  scw{5}=0; scw{6}=0; scw{7}=0;  scw{8}=0;  
            otherwise
        end    
        if scw{sel}==0
            set(sw,'String','');
            set(sd,'String','');
            set(st,'Value',1);
            errordlg('The scale you selected is lower than chosen temporal scale. If you want to use this template you have to downgrade the temporal scale','Invalid scale');
            return;
        end
                    

        switch sel
            case 2 % Morning/Afternoon/Evening/Night (Simple)
               set(sw,'String',num2str([1 0.50 0.75 0.25]));
               set(sd,'String',num2str([7 5 4 8]*scw{2}));
            case 3 % Morning/Afternoon/Evening/Night (Advanced)
               set(sw,'String',num2str([0.8 0.6 0.4 0.5 0.7 0.3 0.1]));
               set(sd,'String',num2str([4 3 3 2 2 2 8]*scw{3}));              
            case 4 % Hour
               set(sw,'String', num2str( (1:24)/24 ,1) ) ;
               set(sd,'String',num2str(ones(1,24)*scw{4}) );   
            case 9 % Year
               if scle==4
                    set(sw,'String', num2str( (1:52)/52 ,1) ) ;
                    set(sd,'String',num2str(floor(ones(1,52))) );   
               elseif scle==5
                    set(sw,'String', num2str( (1:12)/12 ,1) ) ;
                    set(sd,'String',num2str(floor(ones(1,12))) );   
               elseif scle==6
                    set(sw,'String', num2str( (1:4)/4 ,1) ) ;
                    set(sd,'String',num2str(floor(ones(1,4))) );   
               else    
                    set(sw,'String', num2str( (1:365)/365 ,1) ) ;
                    set(sd,'String',num2str(floor(ones(1,365)*scw{5})) );   
               end 
            case 5 % Week
               set(sw,'String', num2str( (1:7)/7 ,1) ) ;
               set(sd,'String',num2str(ones(1,7)*scw{6}) );   
            case 6 % Weekend/Weekends  
               set(sw,'String',num2str([0.5 1]));
               set(sd,'String',num2str([2 5]*scw{7}));
            case 7 % Month
                set(sw,'String', num2str( (1:12)/12 ,1) ) ;
                set(sd,'String', num2str(floor([31 28 31 30 31 30 31 31 30 31 30 31]*scw{8})) );   
            case 8 % Season
               set(sw,'String',num2str([0.75 1 0.5 0.25]));
               set(sd,'String',num2str(floor([91 93 93 88]*scw{9})));            
            otherwise
               set(sw,'String','');
               set(sd,'String','');
        end
        
        status_update();
end

function seasonal_reset(PushButton, EventData)
 	tobj=tab{findtab('SeasonalEffects')};
    ss=findobj(tobj,'Tag','seasonal_scale');
    tbl=findobj(tobj,'Tag','uitable_seasonal');
    set(ss,'Enable','on');
    set(tbl,'Data',{blanks(0),blanks(0),blanks(0)});
    for i=1:1
        RowsN(i,1)={['Factor ',sprintf('%i',i)]};
    end    
    set(tbl,'RowName',RowsN);
    i_defiend_seasonal=0;
    status_update();
end

function seasonal_add(PushButton, EventData)

    R_seasonal_max=i_R_Seasonal;

    if i_seasonal==1 && i_periodic==1 && i_mode_Periodic==i_mode_Seasonal
        if i_defiend_periodic==0
            rst=1;
        else
            rst=i_defiend_periodic;
        end    
        R_seasonal_max=i_R_Seasonal-rst;
    end
	if i_seasonal==1 && i_periodic==1 && i_mode_Periodic~=i_mode_Seasonal
        R_seasonal_max=i_R_Seasonal;
    end
    if i_defiend_seasonal+1>R_seasonal_max
        errordlg(sprintf('The maximum number of factors for seasonal effects should not exceed %i',R_seasonal_max));
        return;
    end
    
	tobj=tab{findtab('SeasonalEffects')};
	sw=findobj(tobj,'Tag','seasonal_weights');
	sd=findobj(tobj,'Tag','seasonal_durations');
    grw=findobj(tobj,'Tag','seasonal_growth');
    tbl=findobj(tobj,'Tag','uitable_seasonal');
    ss=findobj(tobj,'Tag','seasonal_scale');
    dat=get(tbl,'Data');
   
    c1=str2num(get(sw,'String'));
    c2=str2num(get(sd,'String'));
    if length(c1)<2
         errordlg('The length of array should be at least 2');
         return;
    end
    if length(c1)~=length(c2)
         errordlg(sprintf('The length of weight array is %i while length of duration is %i! The size of these two array should be the same.',length(c1),length(c2)));
         return;
    end
    for i=1:length(c2)
        if ~rem(c2(i),1)==0
            errordlg(sprintf('in the duration array only integers are allowed. We found %f as non-integer value! Please revise it. ',c2(i)));
            return;
        end
    end
    
    if all(cellfun(@isempty, dat(:)))    
       dat=[];
       set(tbl,'Data',dat);
    end    
    row{1,1}=get(sw,'String');
    row{1,2}=get(sd,'String');
    row{1,3}=get(grw,'String');
    
    dat=[dat; row];
    set(tbl,'Data',dat);
    set(tbl,'ColumnEditable',logical(zeros(1,size(dat,1))) );
    for i=1:size(dat,1)
        RowsN(i,1)={['Factor ',sprintf('%i',i)]};
    end    
    set(tbl,'RowName',RowsN);
    set(ss,'Enable','off');
    
    i_defiend_seasonal=size(dat,1);
    status_update();
end
    
function periodic_generate(PushButton, EventData)
    tobj=tab{findtab('PeriodWave')};
    
    frequency_from=str2num(get(findobj(tobj,'Tag','frequency_from'),'String'));
    frequency_to=str2num(get(findobj(tobj,'Tag','frequency_to'),'String'));
    nwaves_from=str2num(get(findobj(tobj,'Tag','nwaves_from'),'String'));
    nwaves_to=str2num(get(findobj(tobj,'Tag','nwaves_to'),'String'));
    waveypes={'cos','sin','sq','saw'};
    periodic_typ(1)=get(findobj(tobj,'Tag','periodic_cos'),'Value');
    periodic_typ(2)=get(findobj(tobj,'Tag','periodic_sin'),'Value');
    periodic_typ(3)=get(findobj(tobj,'Tag','periodic_sq'),'Value');
    periodic_typ(4)=get(findobj(tobj,'Tag','periodic_saw'),'Value');
    R_periodic=get(findobj(tobj,'Tag','R_periodic'),'Value');
    waveypes(find(periodic_typ==0))=[];
    waves=waveypes(randi([1 length(waveypes)],1,R_periodic))';
    freqs=randi([frequency_from frequency_to],1,R_periodic)';
    nwaves=randi([nwaves_from nwaves_to],1,R_periodic)';
   
    dat(:,1)=waves;
    for i=1:length(nwaves)
        dat(i,2)={nwaves(i)};
    end
    for i=1:length(freqs)
        dat(i,3)={freqs(i)};
    end
    
   % period_table2=get(findobj(tobj,'Tag','uitable_periodic_table'),'Data')
    set(findobj(tobj,'Tag','uitable_periodic_table'),'Data',dat);
    set(findobj(tobj,'Tag','uitable_periodic_table'),'ColumnEditable',logical(ones(1,R_periodic*3)));

    for i=1:R_periodic
        RowsN(i,1)={['Factor ',sprintf('%i',i)]};
    end    
    set(findobj(tobj,'Tag','uitable_periodic_table'),'RowName',RowsN);
    i_defiend_periodic=length(nwaves);
    status_update();
    
end    
function periodic_edit(PushButton, EventData)
    dat=get(PushButton,'Data');
    for i=1:size(dat,1)
       allfilled(i)=isempty(find([strcmp(dat(i,1),'') strcmp(dat(i,2),'') strcmp(dat(i,3),'')]==1));
    end
    i_defiend_periodic=length(find(allfilled==1));
    status_update();
end
    
function periodic_empty(PushButton, EventData)
    tobj=tab{findtab('PeriodWave')}; 
    R_periodic=get(findobj(tobj,'Tag','R_periodic'),'Value');
    emp{1}='';
    dat(:,1)=emp(ones(R_periodic,1));
    dat(:,2)=dat(:,1);
    dat(:,3)=dat(:,1);
    set(findobj(tobj,'Tag','uitable_periodic_table'),'Data',dat);
    set(findobj(tobj,'Tag','uitable_periodic_table'),'ColumnEditable',logical(ones(1,R_periodic*3)));
    
   for i=1:R_periodic
        RowsN(1,i)={['Factor ',sprintf('%i',i)]};
    end    
    set(findobj(tobj,'Tag','uitable_periodic_table'),'RowName',RowsN);
    i_defiend_periodic=0;
    status_update();
end

function streaming_update(PushButton, EventData)
        tobj=tab{findtab('Streaming')};
        c_streaming=get(findobj(tobj,'Tag','c_streaming'),'Value');
        estrmobj=findobj(tobj,'Tag','e_streaming_param');
        if c_streaming==1
            set(estrmobj,'Enable','on');
            i_streaming=str2num(get(estrmobj,'String'));
        else
            set(estrmobj,'Enable','off');
            i_streaming=0;
        end
        status_update();
end        

function settings_update(PushButton, EventData)
    tobj=tab{findtab('Settings')};
    c_angle=get(findobj(tobj,'Tag','c_angle'),'Value');
    c_correlation=get(findobj(tobj,'Tag','c_correlation'),'Value');
    i_normal=get(findobj(tobj,'Tag','c_normal'),'Value');
    i_fixsign=get(findobj(tobj,'Tag','c_sign'),'Value');

    angleobj=findobj(tobj,'Tag','e_angle');
    corrobj=findobj(tobj,'Tag','e_correlation');
    set(angleobj,'Enable',offon(c_angle));
    set(corrobj,'Enable',offon(c_correlation));

    i_s_angle=0;
    i_s_corr=0;
    if c_angle==1
        i_s_angle=str2num(get(angleobj,'String'));
    end
    if c_correlation==1
        i_s_corr=str2num(get(corrobj,'String'));
    end
       
    status_update();
    
end

% ***************************************************
% Inital settings for the tabs
% ***************************************************

function tab_check(PushButton, EventData)
    switch tgroup.SelectedTab.Title
        case 'PeriodWave'
            tobj=tab{findtab('PeriodWave')};
            rpobj=findobj(tobj,'Tag','R_periodic');
            max_i_R_Periodic=i_R_Periodic;
            if i_defiend_seasonal==0
                rst=1;
            else
                rst=i_defiend_seasonal;
            end    
            if i_seasonal==1 && i_periodic==1 && i_mode_Periodic==i_mode_Seasonal
                max_i_R_Periodic=max_i_R_Periodic-rst;
            end
            if max_i_R_Periodic==0
                max_i_R_Periodic=1;
            end    
            if max_i_R_Periodic>0
                for i=1:max_i_R_Periodic
                    rpobjs(i)={i};
                end    
                set (rpobj,'String',rpobjs);
            end    

            if i_seasonal==0 && i_periodic==1
              set (rpobj,'Value',max_i_R_Periodic);
              set (rpobj,'Enable','off');
            end
            if i_seasonal==1 && i_periodic==1 && i_mode_Periodic==i_mode_Seasonal
              set (rpobj,'Value',ceil(max_i_R_Periodic/2));
              set (rpobj,'Enable','on');
            end
            if i_seasonal==1 && i_periodic==1 && i_mode_Periodic~=i_mode_Seasonal
              set (rpobj,'Value',max_i_R_Periodic);
              set (rpobj,'Enable','off');
            end
            if i_periodic==0
                errordlg('You have not defined any periodic waves in neither of modes!');
                tgroup.SelectedTab = tab{findtab('Modes')};
            end
            if i_defiend_periodic==0
                periodic_empty('','');
            end    
        case 'SeasonalEffects'
            tobj=tab{findtab('SeasonalEffects')};
            if i_seasonal==0 
                errordlg('You have not defined any mdoe for seasonal effecs!');
                tgroup.SelectedTab = tab{findtab('Modes')};
            end
        case 'ChangePoints'
            tobj=tab{findtab('ChangePoints')};
            if i_temporal==0
                errordlg('You have not defined any temporal dimension');
                tgroup.SelectedTab = tab{findtab('Modes')};
            end   
            
            tobj=tab{findtab('ChangePoints')};
            cp_seasonal=findobj(tobj,'Tag','cp_seasonal');
            uitable_changepoints_s=findobj(tobj,'Tag','uitable_changepoints_s');
            changepoints_s_title=findobj(tobj,'Tag','changepoints_s_title');

            cp_periodic=findobj(tobj,'Tag','cp_periodic');
            uitable_changepoints_p=findobj(tobj,'Tag','uitable_changepoints_p');
            changepoints_p_title=findobj(tobj,'Tag','changepoints_p_title');


            if i_seasonal==0
                set(cp_seasonal,'Enable','off');
                set(cp_seasonal,'Value',0);
                set(uitable_changepoints_s,'Visible','off');
                set(changepoints_s_title,'Visible','off');
            else
                set(cp_seasonal,'Enable','on');
                set(uitable_changepoints_s,'Visible','on');
                set(changepoints_s_title,'Visible','on');
            end
            if i_periodic==0
                set(cp_periodic,'Enable','off');
                set(cp_periodic,'Value',0);
                set(uitable_changepoints_p,'Visible','off');
                set(changepoints_p_title,'Visible','off');
            else
                set(cp_periodic,'Enable','on');
                set(uitable_changepoints_p,'Visible','on');
                set(changepoints_p_title,'Visible','on');
            end
        case 'Anomalies'

                if i_added==0
                    errordlg('You have not defined any mode!');
                    tgroup.SelectedTab = tab{findtab('Modes')};
                else
                    tobj=tab{findtab('Anomalies')};
                    anomalies_slide_prev=findobj(tobj,'Tag','anomalies_slide_prev');
                    anomalies_slide_prev.String=mat2str(ceil(i_I*0.1));
                    anomalies_length=findobj(tobj,'Tag','anomalies_length');
                    set (anomalies_length,'Min',0.01);
                    uitable_anomalies=findobj(tobj,'Tag','uitable_anomalies');
                    set(uitable_anomalies,'ColumnEditable',logical([0 0 0 0 0 0]) );
                end
        case 'Generate'
              if i_added==0
                errordlg('You have not defined any mode!');
                tgroup.SelectedTab = tab{findtab('Modes')};
              end      
    end        
end        

% *********************************************************************************************************
% ***********************************      Generating User Interface    ***********************************
% *********************************************************************************************************


tabs={'Modes','PeriodWave','SeasonalEffects','Streaming','ChangePoints','Anomalies','Noise','NonNegative','Sparsity','Settings','Generate','About'};

f = figure('MenuBar','none','Name','SimTensor','NumberTitle','off','Position',[300 300 790 575],'resize','off');
movegui(f,'center');
tgroup = uitabgroup('Parent', f, 'Position',[-0.001 0.1 1.01 0.90],'SelectionChangedFcn',@tab_check);
 
k_start=1;
k_end=length(tabs);

for k=k_start:k_end
    tab{k} = uitab('Parent', tgroup, 'Title', tabs{k});
    if isdeployed
        tab_m_file=fullfile(ctfroot, 'GUIt', [tabs{k} '_export.txt']);
    else
        tab_m_file=fullfile('GUIt', [tabs{k} '_export.txt']);
        if ~exist(tab_m_file)
            copyfile('template.fig',['GUI/',tabs{k},'.fig']);
            copyfile('template_export.m',['GUI/',tabs{k},'_export.m']);
        end
    end
    if exist(tab_m_file)
        txt = fileread(tab_m_file);
        % groups
        ufields=strsplit(txt,['appdata.lastValidTag = ',char(39),'uibuttongroup']);
        h{1+k*100}=tab{k};
        for i=2:length(ufields)
            ps{i-1}=cutstr(ufields{i},'Position');
            tt{i-1}=cutstr(ufields{i},'Title');
            tg{i-1}=cutstr(ufields{i},'Tag');
            ud{i-1}=cutstr(ufields{i},'UserData');
            if ~isnumeric(ud{i-1})
                ud{i-1}=1;
            end
            h{ud{i-1}+k*100} = uibuttongroup('Parent',tab{k},'Title',tt{i-1},'unit','characters','Position',ps{i-1},'Tag',tg{i-1});
        end
        
         % tables
        tfields=strsplit(txt,['appdata.lastValidTag = ',char(39),'uitable']);
        h{1+k*100}=tab{k};
        for i=2:length(tfields)
            ps1{i-1}=cutstr(tfields{i},'Position');
            tg1{i-1}=cutstr(tfields{i},'Tag');
            ud1{i-1}=cutstr(tfields{i},'UserData');
            cfe{i-1}=cutstr(tfields{i},'CellEditCallback');
            if ~isnumeric(ud{i-1})
                ud1{i-1}=1;
            end
            dat1{i-1}=cutstr(tfields{i},'Data');
            colformat{i-1}=cutstr(tfields{i},'ColumnFormat');
            colwidth{i-1}=cutstr(tfields{i},'ColumnWidth');
            if isempty(colwidth{i-1})
                colwidth{i-1}='auto';
            end
            cols1{i-1}=cutstr(tfields{i},'ColumnName');
            rows1{i-1}=cutstr(tfields{i},'RowName');
            vis1{i-1}=cutstr(tfields{i},'Visible');
            vis1{i-1}='on';
            if isempty(vis1{i-1})
                vis1{i-1}='on';
            end

            ht{ud1{i-1}+k*100} = uitable('Parent',h{ud1{i-1}+k*100},'unit','characters','Visible',vis1{i-1},'Position',ps1{i-1},'Tag',tg1{i-1},'Data',dat1{i-1},'ColumnName',cols1{i-1},'ColumnWidth',colwidth{i-1},'ColumnFormat',colformat{i-1},'RowName',rows1{i-1},'ColumnEditable',true(1,length(cols1{i-1})));
            if  cfe{i-1}
              eval(['set(ht{',sprintf('%i',ud1{i-1}+k*100),'}, ',char(39),'CellEditCallback',char(39),',{@', cfe{i-1},'})']);
            end 
        end
        
        
        
        %uicontrol
        s1=strfind(txt,'uicontrol');
        s2=strfind(txt,'function local_CreateFcn');
        selected=txt(s1:s2-1);
        fields=strsplit(selected,'uicontrol(');
        
        
        for i=2:length(fields)
            ttle{i-1}=cutstr(fields{i},'String');
            par{i-1}=cutstr(fields{i},'Parent',1);
            pos{i-1}=cutstr(fields{i},'Position');
            style{i-1}=cutstr(fields{i},'Style');
            if strcmp(style{i-1},'checkbox') || strcmp(style{i-1},'Checkbox') || strcmp(style{i-1},'radiobutton') || strcmp(style{i-1},'popupmenu') || strcmp(style{i-1},'slider')
                val{i-1}=str2num(cutstr(fields{i},'Value',1)); 
                if isempty(val{i-1})
                     val{i-1}=0;
                end     
            end
            fontsize{i-1}=cutstr(fields{i},'FontSize');
            cf{i-1}=cutstr(fields{i},'Callback');
            enbl{i-1}=cutstr(fields{i},'Enable');
            if isempty(enbl{i-1})
                enbl{i-1}='on';
            end
            visib{i-1}=cutstr(fields{i},'Visible');
            if isempty(visib{i-1})
                visib{i-1}='on';
            end
            tags{i-1}= cutstr(fields{i},'Tag');
            usd{i-1}=cutstr(fields{i},'UserData');
            if ~isnumeric(usd{i-1})
                usd{i-1}=1;
            end
            if strcmp(style{i-1},'checkbox') || strcmp(style{i-1},'radiobutton') || strcmp(style{i-1},'popupmenu') || strcmp(style{i-1},'slider')
                fd{i-1+k*100}= uicontrol('Parent', h{usd{i-1}+k*100} , 'HorizontalAlignment','left', 'Value',val{i-1},'Style',style{i-1},'unit','characters','Position',pos{i-1},'string',ttle{i-1},'Tag',tags{i-1},'Enable',enbl{i-1},'Visible',visib{i-1});                
            elseif strcmp(tags{i-1},'edit_multiline')
                fd{i-1+k*100}= uicontrol('Parent', h{usd{i-1}+k*100} ,  'HorizontalAlignment','left', 'Style',style{i-1},'unit','characters','Position',pos{i-1},'string',ttle{i-1},'Tag',tags{i-1},'Enable',enbl{i-1},'Max',10,'Min',1);                
            else    
                fd{i-1+k*100}= uicontrol('Parent', h{usd{i-1}+k*100} ,  'HorizontalAlignment','left', 'Style',style{i-1},'unit','characters','Position',pos{i-1},'string',ttle{i-1},'Tag',tags{i-1},'Enable',enbl{i-1});                
            end    
            if  cf{i-1}
                eval(['set(fd{',sprintf('%i',i-1+k*100),'}, ',char(39),'Callback',char(39),',{@', cf{i-1},'})']);
            end 
            if  fontsize{i-1}
                eval(['set(fd{',sprintf('%i',i-1+k*100),'}, ',char(39),'FontSize',char(39),',',sprintf('%i',fontsize{i-1}),')']);
            end    
            
            
            

        end    
    end   
end


sb = uipanel(...
'Parent',f,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units','characters',...
'Title',{  '' },...
'Position',[-0.02 -0.02 158.5 4.5],...
'ParentMode','manual',...
'Tag','uipanel1');

status_lbl={'Strcutrue','Size','Rank','Tensor Order','Fixed Angle','Fixed correlation','Temporal', 'Periodic',...
'Seasonal','Streaming','Changes','Anomalies','Noise','Non-negative','Sparse','Normalize','Sign Fixing','Category: Regular'};
col_sapce=[0 7 7 7 5 0];
lbl_enbl=zeros(1,length(status_lbl));
lbl_enbl(1:4)=1;
lbl_enbl(end)=1;

hgh=[2.8 1.6  0.4];
for i=1:length(status_lbl)
    hi=mod(i,3);
    if hi==0 
        hi=3;
    end    
    ph=hgh(hi);
    wi=ceil(i/3);
    pw=2.2+(wi-1)*26+col_sapce(wi); 
    prt=strtrim(strsplit(status_lbl{i},':'));
    uicontrol('Parent', sb , 'Style','text','unit','characters','Position',[pw ph 30 1.10],'string',status_lbl{i},'Tag',prt{1},'Enable',offon(lbl_enbl(i)),'HorizontalAlignment','left');                
end    


 
end

