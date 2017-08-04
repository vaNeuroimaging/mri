function [weights1, out_weights, output_save_train, test_output_save_train, AUC,MSE,Annealing_Log] = MLP_PCA_variable_learning_MSE_differential_annealing_v6(training_in, test_in, desired_test, desired_train, weight_ini, num_hid_nodes, num_out_nodes, n_in, num_iter,save_times,stop_point,K,minfun,optnet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function MLP_PCA_variable_learning_MSE_differential_annealing_v6
%
%creates a multilayer perceptron to identify rs-fMRI networks
%
%inputs are of this form:
% training_in & test_in: size = [(#patients*#ROIs) , #voxels ]
% desired_trai & desired_test: size = [ 1 , (#patients*#ROIs) ]
% weight_ini: initial weight for weights1 and out_weights
% num_hid_nodes: number of nodes in the hidden layer
% num_out_nodes: number of nodes in the output layer
% n: learning rate for backward propagation
% num_iter: number of iterations of training
% optnet: 
%
%backpropagation training
%Nicholas Szrama
%10/14/2011
%added variable learning rate trajectories for vector n parameter
%Carl Hacker
%5/28/2012
%added simulated annealing and detailed controls for global optimization 
%Carl Hacker
%8/12/2012
%edited for distribution
%Carl Hacker
%9/25/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Weights
%construct the input matrix as a concatentation of all of the training data
input = training_in;
input = [ones(size(training_in,1),1) , input]; %append a bias term to all inputs
input = double(input);
%initialize weight matrix
weights1 = -weight_ini + 2*weight_ini*rand(size(input,2),num_hid_nodes); %(18612 x num_nodes)
%weights1(:,1) = 0; %make sure the bias nodes do not receive input
out_weights = -weight_ini + 2*weight_ini*rand(num_hid_nodes,num_out_nodes); %one output node ([0,9] describing which class it belongs to)
tic
count=1;
gamma=.2;
SSEpre=inf;AUCpre=-inf;

%% Create Output Matrices
%construct the desired output matrix as a binary matrix (147 by 7)
d_output = double(zeros(size(training_in,1), length(unique(desired_train))));
for ind=1:length(unique(desired_train))
    d_output(desired_train == ind, ind) = 1;
end

d_output_test = double(zeros(size(test_in,1), length(unique(desired_test))));
for ind=1:length(unique(desired_test))
    d_output_test(desired_test == ind, ind) = 1;
end

%% Annealing parameters
out_weights_min=out_weights;
weights1_min=weights1;
Umin=inf;Emin=inf;testAUC_max=-inf;testAUC_max=-inf;SSEmin=inf;
gamma=.12;
max_counter=100;counter=max_counter;release_breakout=0;
Uprevious=inf;Upre=inf;Epre=inf;testAUC_pre=-inf;
T1=K(3);%temperature start
T=T1;
K1=K(1);%nbreakouts
K2=K(2);%attempts
breakthresh=K(4);
breakthreshpos=K(5);
coolrate=K(6);
A=K(7); %the std of noise relative to weights
k1=0;k2=0;
WeightEnergy=1;
hit_count=1;
Annealing_Log=[];
%% define logistic activating function parameters (hyperbolic tangent function)
a = 1;
b = 1/10;
%% construct test data
test_data = test_in;
test_data = [ones(size(test_in,1),1) , test_data]; %append a bias term to all inputs
test_data = double(test_data);

%% save output matrix for training set
% output_save_train = double(zeros(size(input,1),num_out_nodes,num_iter/save_interval+1));
output_save_train = double(zeros(size(input,1),num_out_nodes,length(save_times)));
% test_output_save_train = double(zeros(size(test_data,1),num_out_nodes,num_iter/save_interval+1));
test_output_save_train = double(zeros(size(test_data,1),num_out_nodes,length(save_times)));

%% construct AUC matrix
% AUC = zeros(sum(mod(1:num_iter,save_interval)==0)+1, num_out_nodes);

%% Compute learning rate vector
%unless given as parameter
t1=1;t2=num_iter;
tinds=t1:t2;
if numel(n_in)==1
    nvec=0;
elseif numel(n_in)==3
    x1=n_in(1);x2=n_in(2);k=n_in(3);
    a1=(x2-x1)/(exp(t2/k)-exp(t1/k));
    c1=x1-a1*exp(t1/k);
    nvec=a1*exp(tinds./k)+c1;
    %     figure,plot(nvec)
elseif numel(n_in)==6
    vvo=n_in(1);
    QQo=n_in(2);
    AAo=n_in(3);
    KKo=n_in(4);
    BBo=n_in(5);
    MMo=n_in(6);
    nvec=n_in;
    xCDF=linspace(-1,1,length(save_times));
    nvec=AAo+(KKo-AAo)./(1+QQo*exp(-BBo*(xCDF-MMo))).^(1/vvo);
    nvec=interp1(save_times,nvec,1:num_iter);
    %     figure(1),clf,plot(tinds,nvec),set(gca,'fontsize',12,'xscale','log'),xlim([10,max(tinds)]),pause
%final learning rate vector was an input
elseif numel(n_in)==num_iter
    nvec=n_in;
else
    error('nvec wrong length')
end
%% Perceptron Algorithm
iter=0;
while iter < num_iter
    iter=iter+1;
    %forward propagation
    if numel(nvec)>1
        n=nvec(iter);
    else
        n=n_in;
    end
    %input to hidden layer 1
    V1 = input*weights1; %(147 x num_hid_nodes matrix): each row is a double trial, each column represents the V(j)
    y = a*tanh(b*V1); %y values for hidden layer
    %     y = 1./(1 + exp(-a*V1));
    y(:,1) = 1; %apply bias term
    

    %hidden layer 1 to output
    V3 = y*out_weights; %(147 x num_out_nodes matrix): each row is a double trial, each column represents the V(j)
    %     output = a*tanh(b*V3); %output values (147 by num_out_nodes matrix)
    output = 1./(1 + exp(-a*V3));
    
    %error calculation
    e  = d_output - output; %(147 by num_out_nodes matrix)
    
    %local gradient for output node (d3) and hidden layer 2 (d2) and hidden
    %layer 1 (d1)
    %     phi_prime3 = b/a*(a - output).*(a + output);
    phi_prime3 = a*output.*(1 - output);
    d3 = e.*phi_prime3; %(147 by num_out_nodes matrix)
    
    
    phi_prime1 = b/a*(a - y).*(a + y);
    %     phi_prime1 = a*y.*(1 - y);
    d1 = d3*out_weights'.*phi_prime1; %(147 x num_hid_nodes matrix)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %backward propagation
    
    %update weights to output layer
    out_weights = out_weights + n*y'*d3;
    
    %update weights to hidden layer
    weights1 = weights1 + n*input'*d1;
    
    %ensure that the bias terms weights remain at 1
    y(:,1) = 1;
    
    %perform ROC-AUC analysis at save_times
    if any(iter==save_times) %((mod(iter,save_interval) == 0) || (iter==1))
        
        %save output every save_interval iterations
        %         output_save_train(:,:,floor(iter/save_interval)+1) = output;
        output_save_train(:,:,find(iter==save_times)) = output;
        %forward propagation
        %input to hidden layer 1
        test_V1 = test_data*weights1; %(10K x num_hid_nodes matrix): each row is a double trial, each column represents the V(j)
        test_y = a*tanh(b*test_V1); %y values for hidden layer
        test_y(:,1) = 1; %apply bias term
        % test_y = 1./(1 + exp(-a*test_V1));
        
        %hidden layer 1 to output
        test_V3 = test_y*out_weights; %(10K x num_out_nodes matrix): each row is a double trial, each column represents the V(j)
        %     output = a*tanh(b*V3); %output values (10K by num_out_nodes matrix)
        test_output = 1./(1 + exp(-a*test_V3));
        %         test_output_save_train(:,:,floor(iter/save_interval)+1) = test_output;
        test_output_save_train(:,:,find(iter==save_times)) = test_output;
        
        
        %         MSE_train(find(iter==save_times))=norm(e(:));
        e  = d_output - output;
        e_test  = d_output_test - test_output;
        %         MSE_test(find(iter==save_times))=norm(e_test(:));
        nnets=7;
        for et = 1:nnets+2
            
            if et < nnets+1
                etn=e(:,et);
                etn_test=e_test(:,et);
            end
            if et == nnets+1
                etn=e;
                etn_test=e_test;
            end
            if et == nnets+2
                etn=e(:,1:7);
                etn_test=e_test(:,1:7);
            end
            MSE{et}(find(iter==save_times),:)=[norm(etn(:))/sqrt(numel(etn(:))),norm(etn_test(:))/sqrt(numel(etn_test(:)))];
        end
        %determine thresholding values
        TPrate = zeros(101, num_out_nodes);
        FPrate = zeros(101, num_out_nodes);
        ROC_cdf_y_vals=linspace(0,1,100)';
        for net_ind=1:num_out_nodes
            [n1, x1] = hist(test_output(:, net_ind), 500);
            x1=[0 x1 1];
            netcdf = cumsum(n1)/sum(n1);
            netcdf=[0 netcdf 1];
            
            %find number of cdf x inds with value less than test y index
            %             xroc_defining_indices=sum(netcdf(ones(100,1),:) < ROC_cdf_y_vals(:,ones(size(netcdf,2),1)),2);
            %             xroc=x1(xroc_defining_indices);
            xroc_defining_indices=sum(netcdf(ones(100,1),:) < ROC_cdf_y_vals(:,ones(size(netcdf,2),1)),2);
            xroc_defining_indices=xroc_defining_indices(xroc_defining_indices>0);
            xroc=x1([ 1;xroc_defining_indices;length(x1) ] );
            
            for xind=1:length(xroc)
                TPrate(xind, net_ind) = sum((test_output(:, net_ind) >= xroc(xind)) & (desired_test == net_ind)')/(sum((test_output(:, net_ind) >= xroc(xind)) & (desired_test == net_ind)') + sum((test_output(:, net_ind) < xroc(xind)) & (desired_test == net_ind)' ));
                FPrate(xind, net_ind) = sum((test_output(:, net_ind) >= xroc(xind)) & (desired_test ~= net_ind)')/(sum((test_output(:, net_ind) >= xroc(xind)) & (desired_test ~= net_ind)') + sum((test_output(:, net_ind) < xroc(xind)) & (desired_test ~= net_ind)'));
            end
            AUC(count, net_ind) = -trapz(FPrate(:,net_ind), TPrate(:,net_ind));
        end
        
        %Annealing breakout code
        %Definition of objective functions  (AUC or SSE)
        if optnet==0
            esource=e_test(:,1:7);
            testAUC=mean(AUC(count, 1:7),2);
        else
            esource=e_test(:,optnet);
            testAUC=AUC(count, optnet);
        end
        SSE=sum(esource(:).^2);
        RMS=sqrt(SSE/numel(esource));
        E_training= e(:,1:7); %#ok<*NASGU>
        Etr=sum(E_training(:).^2);
        RMStr=sqrt(Etr/numel(E_training));
        U=SSE+gamma*(sum(weights1(:).^2)+sum(out_weights(:).^2));
        Uweights=(sum(weights1(:).^2)+sum(out_weights(:).^2));
        rel_error=abs(SSE-SSEpre)/(.5*(SSE+SSEpre+1e-10));
%         testAUC=mean(AUC(count, 1:7),2);
        esource_all=e_test(:,1:7);
        testAUC_all=mean(AUC(count, 1:7),2);
        SSE_all=sum(esource_all(:).^2);
        RMS_all=sqrt(SSE_all/numel(esource_all));
        
        Breakout_Logic=0;
        switch upper(minfun)
            case 'SSE'
                if ~isnan(rel_error)&&(rel_error<breakthresh)||(rel_error<breakthreshpos&&SSE>SSEpre)
                    Breakout_Logic=1;
                end
            case 'AUC'
                if ~isnan(rel_error)&&(rel_error<breakthresh)||(rel_error<breakthreshpos&&testAUC<testAUC_pre&&testAUC>.8)
                    Breakout_Logic=1;
                end
            case 'NONE'
                Breakout_Logic=0;
            otherwise
                error('No objective function given')
        end
        
        if ~Breakout_Logic
            testAUC_pre=testAUC;
            Upre=U;
            SSEpre=SSE;
        end
%         
%         fprintf('%-.2f AUC%2.4f U:%2.1f Uw:%2.1f Etr:%3.1f(%2.3f) SSE:%3.1f (%2.5f) k1=%d/%d k2=%d/%d opt:%d%s\n',...
%             iter/num_iter,mean(AUC(count, 1:7)),Uweights+Etr,Uweights,Etr,sqrt(Etr)./sqrt(numel(e)),SSE,sqrt(SSE)./sqrt(numel(esource)),k1,K1,k2,K2,optnet,minfun)
%         
        fprintf('%d AUC_U:%2.4f AUC:%2.4f Uw:%2.1f Etr:%3.2f(%2.5f) E:%3.1f(%2.5f) U:%3.1f(%2.5f) k1=%d/%d k2=%d/%d opt:%d%s\n',...
    fix(iter/num_iter*100),testAUC,testAUC_all,Uweights,Etr,RMStr,SSE_all,RMS_all,SSE,RMS,k1,K1,k2,K2,optnet,minfun)

    if all(Breakout_Logic)
            %% Attempt breakout
            if k1<K1
                toc
%                 iter=1000;
                k1=k1+1;Uprevious=inf;counter=max_counter;%6;release_breakout=0;
                iter=10; %reset counter during bump
                A_Hit=0;
                switch upper(minfun)
                    case 'SSE'
                        if SSE<SSEmin
                            A_Hit=1;
                        end
                    case 'AUC'
                        if testAUC>testAUC_max
                            A_Hit=1;
                        end
                    otherwise
                        error('No objective function given')
                end
                if A_Hit
                    fprintf('change old SSEmin %f new %f\n',SSEmin,SSE)
                    SSEmin=SSE;testAUC_max=testAUC;
                    weights1_min=weights1;
                    out_weights_min=out_weights;
                    Annealing_Log(hit_count,:)=[testAUC,testAUC_all,Uweights,Etr,SSE_all,RMS_all,SSE,RMS];
                    hit_count=hit_count+1;
                end
                
                % Cooling function
                T=A.*(coolrate^(k1+k2));
                % Hyperbolic weight perturbation
                WeightEnergy=1;%Mean change in weights
%                 y1=double(solve(sprintf('(2*log(%f + 1) - %f)-(2*log(yy + 1) - yy)=%f*(%f-yy)',T,T,WeightEnergy,T)));
                y1=fzero(@(a)1/(T-a)*(T - a - 4*log(T + 1) + 4*log(a + 1) - 4/(T + 1) + 4/(a + 1))-WeightEnergy,-.2);

                tlow=y1;thigh=T;
                Urand=tlow+(thigh-tlow)*rand(size(weights1_min));range(Urand(:));
                weights1=weights1_min.*(1-Urand)./(1+Urand);
                Urand=tlow+(thigh-tlow)*(rand(size(out_weights_min)));
                out_weights=out_weights_min.*(1-Urand)./(1+Urand);
                fprintf('\t \t bump weights:\n\tcurrent AUCmax: %3.4f SSEmin: %3.1f \n\tT:%2.3f\n',testAUC_max,SSEmin,T)
            else
                if k2<K2
                    k1=1;
                    k2=k2+1;
                else
                    break;
                end
            end
            tic
            
            
        end
        %         fprintf('%-.2f%%-Trn:%2.1f Tst:%2.1f Tst1-7:%2.4f AUC%2.4f delE%2.4f Uw:%2.1f SSE:%3.1f\tre:%2.5f k1=%d-K1\n',...
        %             iter/num_iter,MSE{nnets+1}(end,1),MSE{nnets+1}(end,2),MSE{nnets+2}(end,2),mean(AUC(count, 1:7)),mean(diff(MSE{nnets+2}(max(end-5,1):end,2)),1),Uweights,SSE,rel_error,k1)
        pause(.00000001)
        count = count + 1;
    end
    if stop_point>0
        if iter==stop_point
            break
        end
    end
    %     0.197
end

%save final output
output_save_train(:,:,end) = output;



