% addpath('../common');
% Parameters Setting
param = SetSCDLParams();
nClass = 32;
param.win = 5;
param.K = 256;
lassoParam=param.lassoParam;
% x01=imread('20041126_TM.tif');
% y01=imread('MOD09GA_A20041126.tif');
% x02=imread('20041212_TM.tif');
% y02=imread('MOD09GA_A20041212.tif');
% A01=imcrop(x01,[900 1100 999 999]);
%  A02=imcrop(y01,[900 1100 999 999]);
% A04=imcrop(x02,[900 1100 999 999]);
%  A03=imcrop(y02,[900 1100 999 999]);
%  x01=imread('L7SR.05-24-01.r-g-nir.tif');
% y01=imread('MOD09GHK.05-24-01.r-g-nir.tif');
% x02=imread('L7SR.07-11-01.r-g-nir.tif');
% y02=imread('MOD09GHK.07-11-01.r-g-nir.tif');
%    A01=imcrop(x01,[1 200 999 999]);
%  A02=imcrop(y01,[1 200 999 999]);
%  A04=imcrop(x02,[1 200 999 999]);
%  A03=imcrop(y02,[1 200 999 999]);
% A01=imread('lan2\L_20020111.tif');
%  A02=imread('MOD09GA_A20020112.tif');
%  A04=imread('L_20020212.tif');
% A03=imread('MOD09GA_A20020213.tif');
% x01=imread('Landsat20041126.tif');
% y01=imread('MODIS20041126.tif');
% x02=imread('Landsat20041228.tif');
% y02=imread('MODIS20041228.tif');
%  psf = fspecial('gauss',16, 8);
%  for i=1:3
% A02(:,:,i)=conv2(double(A02(:,:,i)),psf,'same');
% A03(:,:,i)=conv2(double(A03(:,:,i)),psf,'same');
%  end
A01=imread('L020125.tif');
A02=imread('M020125.tif');
A03=imread('M020226.tif');
A04=imread('trueL026.tif');
A02=imresize(imresize(A02,0.2,'bicubic'),5,'bicubic');
A03=imresize(imresize(A03,0.2,'bicubic'),5,'bicubic');
A1=double(A01)-double(A02);
A2=double(A02)-imresize(imresize(double(A02),0.2,'bicubic'),5,'bicubic');
A3=double(A03)-imresize(imresize(double(A03),0.2,'bicubic'),5,'bicubic');
% A011=imresize(double(A01),0.25,'bicubic');
% A022=imresize(double(A02),0.25,'bicubic');
% A033=imresize(double(A03),0.25,'bicubic');
% A1=A011-A022;
% A2=A022-imresize(imresize(A022,0.25,'bicubic'),4,'bicubic');
% A3=A033-imresize(imresize(A033,0.25,'bicubic'),4,'bicubic');
[m,n]=size(A1(:,:,1));
k=0;
patch_size=param.win;
c1=m/patch_size;
for i11=1:3
for i = 1:c1
for j = 1:c1
m_start=1+(i-1)*fix(m/c1);
m_end=i*fix(m/c1);
n_start=1+(j-1)*fix(n/c1);
n_end=j*fix(n/c1);
A1A=A1(m_start:m_end,n_start:n_end,:);%将每块读入矩阵
k=k+1;
AA1=A1A(:,:,i11);
A2A=A2(m_start:m_end,n_start:n_end,:);%将每块读入矩阵
AA2=A2A(:,:,i11);
A3A=A3(m_start:m_end,n_start:n_end,:);%将每块读入矩阵
AA3=A3A(:,:,i11);
L(:,k)=AA1(:);
L2(:,k)=AA2(:);
L3(:,k)=AA3(:);
end
end
var1=var(double(L(:)));
var2=var(double(L2(:)));
var3=var(double(L3(:)));
mean1=mean(mean(double(L)));
mean2=mean(mean(double(L2)));
mean3=mean(mean(double(L3)));
L = (double(L) - repmat(mean(mean(double(L))), [param.win^2 1]))./(repmat(var1, [param.win^2 1]));
L2 = (double(L2) - repmat(mean(mean(double(L2))), [param.win^2 1]))./(repmat(var2, [param.win^2 1]));
L3 = (double(L3) - repmat(mean(mean(double(L3))), [param.win^2 1]))./(repmat(var3, [param.win^2 1]));
Dl= mexTrainDL(double(L2), param.lassoParam);
Al=mexLasso(L2,Dl,lassoParam);
Dh=L*full(Al)'*pinv(full(Al)*full(Al)');
Al2=mexLasso(L3,Dl, lassoParam);
Xh1=Dh*Al.*repmat(var2, [param.win^2 1]);
Xh2=Dh*Al2.*repmat(var3, [param.win^2 1]);
Xh11=Xh1+mean2;
Xh22=Xh2+mean3;
xh1=Xh11;
xh2=Xh22;
XH1=cal18(xh1,param.win)+imresize(imresize(double(A02(:,:,i11)),0.25,'bicubic'),4,'bicubic');
XH2=cal18(xh2,param.win)+imresize(imresize(double(A03(:,:,i11)),0.25,'bicubic'),4,'bicubic');
XH(:,:,i11)=(XH2+60)./(XH1+60).*double(A01(:,:,i11));
k=0;

% L4=zeros(2*patch_size*patch_size,(m/patch_size)*(m/patch_size));
% 
% var1=var(double(L(:)));
% var2=var(double(L2(:)));
% var3=var(double(L3(:)));
% XH_t = (double(L) - repmat(mean(mean(double(L))), [param.win^2 1]))./(repmat(var1, [param.win^2 1]));
% 	XL_t = (double(L2) - repmat(mean(mean(double(L2))), [param.win^2 1]))./repmat(var2, [param.win^2 1]);
% L4(1:patch_size*patch_size,:)=XH_t;
% L4(patch_size*patch_size+1:2*patch_size*patch_size,:)=XL_t;
% D = mexTrainDL(double(L4), param.lassoParam);
% Dh = D(1:param.win^2,:);
% 	Dl = D(param.win^2+1:end,:);
%     xk=repmat(mean(mean(double(L3))), [param.win^2 1]);
%     X_t=(double(L3) - repmat(mean(double(L3)), [param.win^2 1]))./repmat(var3, [param.win^2 1]);
%      alphaH = mexLasso(double(XH_t), Dh, lassoParam);
%      kui=double(eye(size(Dh,2)));
%        alphaL = mexLasso([double(X_t);4*full(alphaH)], [Dl;4*kui], lassoParam);
% % alphaL = mexLasso(double(X_t), Dl, lassoParam);
%      xh=Dh*alphaL.*repmat(var3, [param.win^2 1]);
%      xh1=1.5*xh+xk;
% Xh=cal18(xh1,patch_size);
% k=0;
% XH(:,:,i11)=Xh;
end
figure(4)
subplot(133)
imshow(uint8(XH),[])
title('融合结果图');
subplot(131)
imshow(uint8(A01),[])
title('前一时刻图像');
subplot(132)
imshow(A04,[])
title('真实图像');
