function C = TensorTimes(A,B)
% 快速计算 C（：，：，i） = A(:,:,i)*B(:,:,i);
% 大量同大小的矩阵同时乘积，利用张量对每层做矩阵乘积，先把张量升到四维再做点积再sum回三维
A = permute(A,[2,4,3,1]);
C = bsxfun(@times, A,B) ;
C = permute(sum(C,1),[4,2,3,1]);
end