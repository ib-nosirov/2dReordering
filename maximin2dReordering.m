% Author: Ibrohim Nosirov
% Date: 2023-01-15
% Version: 1.0
% Description: a maximin algorithm for ordering points in 2 dimensions into a
% HODLR-type matrix. The algorithm is based on the idea that faraway
% interactions have near-zero, or even slightly negative values, and can
% therefore be compressed into HODLR-type if these interactions reside on the
% off-diagonal blocks of the matrix.
% Input: a set of points in 2 dimensions
% Output: a HODLR-type matrix

% Generate set of 2D points as an array of size 2xN.
% The points are generated randomly in the unit square.
clear;
close all;
N = 100;
dim2Points = [randn(2,100);ones(1,100)];
% Create a radial-basis function kernel
rbfKernel = @(x,y) exp(-norm(x-y)^2);
% Evaluate the kernel at all pairs of points
K = zeros(N,N);
HODLR_Mtrx = K;
for ii = 1:N
	for jj = 1:N
		K(ii,jj) = rbfKernel(dim2Points(1:2,ii),dim2Points(1:2,jj));
	end
end

% Plot the kernel matrix
%figure(1)
%imagesc(K)

validPoints = dim2Points(1:2,dim2Points(3,:) == 1);
validPtIndices = find(dim2Points(3,:) == 1);
focus = validPoints(1:2,end);
for ii = 1:N-1
    % pick the last point from dim2Points and remove it from list of distances.
    dim2Points(3,validPtIndices(end)) = 0; % TODO: fix this index
    % compute the distances of all the other points from the locus
    [~,distIdx] = sort(vecnorm(validPoints(1:2,:) - locus));
    % sort the points by distance
    validPoints = validPoints(:,distIdx);
    % put this vector inside HODLR_Mtrx
    for jj = 1:N-(ii-1) % TODO: this needs to make a lower triangular matrix
		HODLR_Mtrx(jj,ii) = rbfKernel(validPoints(1:2,jj),dim2Points(1:2,ii));
    end
    focus = validPoints(1:2,end);
    validPoints = validPoints(1:2,1:N-ii)
end

figure(2)
imagesc(HODLR_Mtrx)

% Algorithm:
% 