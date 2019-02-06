function [PSNRvector,SSIMvector,avsam1,MQresult] = evaluate(OriData3,output_image,M,N)
p = size(OriData3,3);
PSNRvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);

    I=255*output_image(:,:,i);

      PSNRvector(1,i)=PSNR(J,I,M,N);
end
% dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc');
% 
PSNRvec = mean(PSNRvector);
SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);
%     Jnoise=oriData3_noise(:,:,i);
    I=255*output_image(:,:,i); 
%      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
      [ SSIMvector(1,i),ssim_map] = ssim(J,I);
end
SSIMvec=mean(SSIMvector);
% dlmwrite('SSIMvector.txt',SSIMvector,'delimiter','\t','newline','pc');
sum1=0.0;
for i=1:M;
    for j=1:N;
       T=OriData3(i,j,:);
       T=T(:)';
       H=output_image(i,j,:);
       H=H(:)';
       sum1=sum1+SAM(T,H);
    end
end
avsam1=sum1/(M*N);
MQresult = zeros(p,1);
% for i =1:p;
% AnisoSet=AnisoSetEst(reshape(output_image(:,:,i),M,N),8);
% MQresult(i) = MetricQ(reshape(output_image(:,:,i),M,N),8,AnisoSet);
% end
% MQvec = mean(MQresult);