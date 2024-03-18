function [Q,Sigma] = TakagiSVD(A)

% Description: Perform Takagi's factorization of a complex AND symmetric nxn matrix A such
% that A = Q Sigma Q^T, where Q is unitary and Sigma = diag(sigma1, sigma2,...) with
% sigma1>=sigma2>=....>=0 In fact, Sigma contains the singular values of A
%
% We use an SVD-based algorithm.
% Input parameters:
% A = nxn complex symmetric matrix
%
% Output parameters:
% Q: unitary matrix
% Sigma: diagonal matrix with Takagi's values (singular values)
%
% Ignacio Santamaria, UC March 2024

[n,m] = size(A);
if (n~=m)||(~issymmetric(A))   % if A is non-square or non symmetric
    disp('Error: matrix must be symmetric and square to perform Takagi factorization');
    return
end

if norm(A'*A-eye(n))<1e-10
    Q = sqrtm(A);
    Sigma = eye(n);
else
    [F_SVD,Sigma,G_SVD] = svd(A);
    dummy = F_SVD'*conj(G_SVD);
    phases = angle(diag(dummy))/2;    % Eigenvectors post-correction
    Q = F_SVD*diag(exp(1i*phases));
    % A = Q*Sigma*Q.';
end
