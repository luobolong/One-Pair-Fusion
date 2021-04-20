function cal=cal18(Xl,wi)
% load dictionary 
% Dh=xlsread('D:\Users\0x8A63F77D\Desktop\46\W.xlsx');

Dh=Xl;
patch_size=wi;
nRow=500/wi;
nCol=500/wi;
D=Dh';
w=nCol*patch_size;
h=nRow*patch_size;
gridx = 1:patch_size :w-patch_size;  
gridx = [gridx, w-patch_size+1];  
gridy = 1:patch_size : h-patch_size;  
gridy = [gridy, h-patch_size+1];  
K=nRow*nCol; %字典原子总数
DD=cell(1,K);
row=length(gridx);
col=length(gridy);
hIm=zeros([w,h]);
for i=1:K
    DD{i}=D(i,:);
end
for ii = 1:length(gridx)
    for jj = 1:length(gridy)
        yy = gridx(ii);
        xx = gridy(jj);
        if (ii-1)*nRow+jj >K
            break
        end
        temp=DD{(ii-1)*nCol+jj};
        hPatch=reshape(temp,[patch_size,patch_size]);
        hIm(yy:yy+patch_size-1, xx:xx+patch_size-1) = hIm(yy:yy+patch_size-1, xx:xx+patch_size-1) +hPatch;
     end
end
cal=hIm;
end
%  hIm=medfilt2(hIm,[5 5]);
% 
% hIm=imresize(hIm,4,'bicubic');
% for i=1:400
%     for j=1:400
% if hIm(i,j)<0
%     hIm(i,j)=-hIm(i,j);
% end
% % % if hIm(i,j)<3
% % %     hIm(i,j)=60;
% % % end
% if hIm(i,j)>255
%     hIm(i,j)=255;
% end
%     end
% end
% for i=390:398
%     for j=1:400
% hIm(i,j)=255;
%     end
% end
% figure;
% subplot(221)
% imshow(uint8(hIm),[])
% colormap(gray);
% axis image;
% subplot(222)
% imshow(P1(:,:,1))
% subplot(223)
% imshow(p3(:,:,1))
% subplot(224)
% P4=P4(:,:,1);
% P4=imresize(P4,0.5);
% imshow(P4)
%  P2= imread('D:\Users\0x8A63F77D\Desktop\tu\modis0524.tif');
%  imshow(p12(:,:,1))