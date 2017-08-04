function eigvals_per = PCA_reduction_cov_onecomp(covmat)
%TOL

%Calculate covariance of data
%covmat = cov(mat);

%Calculate eigenvectors and eigenvalues
[V D] = eig(covmat);

eigvals = diag(D);

[eigvals_sort ind] = sort(eigvals,'descend');

eigvals_per = (eigvals_sort./sum(eigvals_sort))*100;

            
            

