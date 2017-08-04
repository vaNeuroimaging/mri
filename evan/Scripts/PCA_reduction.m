function [Y V_use eigvals_per] = PCA_reduction(mat,reductype,reducnum)
%[Y V_use eigvals_per] = PCA_reduction(mat,reductype,reducnum)
%TOL

%Calculate covariance of data
covmat = cov(mat);

%Calculate eigenvectors and eigenvalues
[V D] = eig(covmat);

eigvals = diag(D);

[eigvals_sort ind] = sort(eigvals,'descend');

V_sort = V(:,ind);

eigvals_per = (eigvals_sort./sum(eigvals_sort))*100;

switch reductype
    
    case 'comps'
    
        V_use = V_sort(:,1:reducnum);
        
    case 'percentvar'
        
        count = 1;
        persum = 0;
        while persum < reducnum
            
            persum = persum + eigvals_per(count);
            count = count + 1;
        end
        
        V_use = V_sort(:,1:count);
end

%Dimensionality reduction
Y = V_use'*mat';
            
            

