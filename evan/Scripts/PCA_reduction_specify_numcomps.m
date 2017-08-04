function [Y V_use eigvals_per] = PCA_reduction_specify_numcomps(mat,reducnum)
%[Y V_use eigvals_per] = PCA_reduction(mat,reductype,reducnum)
%TOL

%Calculate covariance of data
covmat = cov(mat);

%Calculate eigenvectors and eigenvalues
[V D] = eigs(covmat,reducnum);

eigvals = diag(D);

[eigvals_sort ind] = sort(eigvals,'descend');

V_use = V(:,ind);

eigvals_per = (eigvals_sort./sum(eigvals_sort))*100;

Y = V_use'*mat';
            
            

