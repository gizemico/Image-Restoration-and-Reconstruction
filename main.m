%% For image-1
clear all;clc;

%read the image
I=imread('noisy1.png');
figure(1);imshow(I);title('noisy1.png');

%take the histogram of image
hst=imcrop(I);
figure(2);imshow(hst);
hist=imhist(hst);
figure(3);bar(hist);

%compute the shifted DFT of the image using functions fft2 and fftshift
F=fftshift(fft2(I));

%create the notchpass and reject filters
[m n]=size(I);
[hnp,hnr]=NotchPassAndRejectFilters('ideal',m,n,120,70,0,2);%d0=120,uk=70,vk=0,N=2.

%filter the image by multiplying the filter with shifted DFT of the image
G=hnp.*F;
R=hnr.*F;

%compute the inverse DFT using ifft2 and abs functions
G2=abs((ifft2((G))));
R2=abs((ifft2((R))));

figure(4);edge(I);title('noisy edge');
figure(5);edge(G2);title('filtered edge');

figure(6);
subplot(1,5,1);imshow(abs(log(F)),[]);title('fourier transform1');
subplot(1,5,2);imshow(hnp,[]);title('notchpass1');
subplot(1,5,3);imshow(abs(log(G)),[]);title('notchpass');
subplot(1,5,4);imshow(abs(log(G)),[]);title('reject filter');
subplot(1,5,5);imshow(G2,[]);title('recovered1');

G2 = uint8(255 * mat2gray(G2));
imwrite(G2,'recovered1.png');



%% For image-2

clear all;clc;

%read the image
I=imread('noisy2.png');
figure(1);imshow(I);title('noisy2.png');

%take the histogram of image
hst=imcrop(I);
figure(2);imshow(hst);
hist=imhist(hst);
figure(3);bar(hist);

%compute the shifted DFT of the image using functions fft2 and fftshift
F=fftshift(fft2(I));

%create the bandpass and reject filters
[m n]=size(I);
[hbp,hbr]=BandPassAndRejectFilters('gaussian',m,n,180,2,180);%d0=180,N=2,w=180.

%filter the image by multiplying the filter with shifted DFT of the image
G=hbp.*F;
R=hbr.*F;

%compute the inverse DFT using ifft2 and abs functions
G2=abs((ifft2((G))));
R2=abs((ifft2((R))));

figure(4);edge(I);title('noisy edge');
figure(5);edge(R2);title('filtered edge');

figure(6);
subplot(1,5,1);imshow(abs(log(F)),[]);title('fourier transform2');
subplot(1,5,2);imshow(hbp,[]);title('bandpass2');
subplot(1,5,3);imshow(abs(log(G)),[]);title('bandpass');
subplot(1,5,4);imshow(abs(log(R)),[]);title('reject filter');
subplot(1,5,5);imshow(R2,[]);title('recovered2');

R2 = uint8(255 * mat2gray(R2));
imwrite(R2,'recovered2.png');

 
    

%% For image-3

clear all;clc;

%read image
I=imread('noisy3.png');
figure(1);imshow(I);title('noisy3.png');

%take the histogram of image
hst=imcrop(I);
figure(2);imshow(hst);
hist=imhist(hst);
figure(3);bar(hist);

%size of image
[m n]=size(I);

%find mean and variance of the block to estimate gaussian
mu=mean(hst(:));
s=std(double(hst(:)));
x=[0:255];hold on;

%form the gaussian
pz=(1/(sqrt(2*pi)*s)) *exp(-((x-mu).^2)/(2*s.^2));

ng=pz*sum(hist)/sum(pz);

plot(x,ng,'r-','linewidth',2);

%padding.
pad=uint8(zeros(m+2,n+2));
pad(2:m+1,2:n+1)=I; %starts with 2, because after enlarge with pad our image will start pixel of (2,2).

%create matrices.
arr=zeros(3,3);
output=uint8(zeros(m,n));
new=zeros(1,9);

for i=2:1:m+1 %starts with 2, and ended with row+1 because after padding, original image inside pad. 
     for j=2:1:n+1 %starts with 2, and ended with column+1 because after padding, original image inside pad. 
         
         summ=0;
         
         arr(1,1)= pad(i-1,j-1);
         arr(2,1)= pad(i,j-1);
         arr(3,1)= pad(i+1,j-1);
         arr(1,2)= pad(i-1,j);
         arr(2,2)= pad(i,j);
         arr(3,2)= pad(i+1,j);
         arr(1,3)= pad(i-1,j+1);
         arr(2,3)= pad(i,j+1);
         arr(3,3)= pad(i+1,j+1);

         x=sort(arr(1,:)); %sort for first row.
         y=sort(arr(2,:)); %sort for second row.
         z=sort(arr(3,:)); %sort for third row.
         new=[x,y,z]; %combine the sorted rows(x,y and z).
         s=sort(new); %sort the new combining array.

         summ=summ+(1/s(1))+(1/s(2))+(1/s(3))+(1/s(4))+(1/s(5))+(1/s(6))+(1/s(7))+(1/s(8))+(1/s(9));%find the mean value.
         harman=(9/summ);%find the harmonic filter.
         output(i-1,j-1)=harman;
         
     end %harmonic filter is usefull for salt and gaussian noise.
end

figure(4);edge(I);title('edge of noisy3');
figure(5);edge(output);title('edge of recovered3');
figure(6);imshow(output);title('recovered3');

imwrite(output,'recovered3.png');


%% For image-4

clear all;clc;

%read image
I=imread('noisy4.png');
figure(1);imshow(I);title('noisy4.png');

%take the histogram of image
hst=imcrop(I);
figure(2);imshow(hst);
hist=imhist(hst);
figure(3);bar(hist);

%size of image
[m n]=size(I);

a=2;b=10;%we assume the values of a and b.
z=randi(255,m,n);% take the value of z as random.

%form the uniform noise
if ((a<=z<=b))
    pz=1/(b-a);
else
    pz=0;
end

ng=pz*sum(hist)/sum(pz);
x=[0:255];hold on;
plot(x,ng,'r-','linewidth',2);

%padding.
pad=uint8(zeros(m+2,n+2));
pad(2:m+1,2:n+1)=I; %starts with 2, because after enlarge with pad our image will start pixel of (2,2).

%create matrices.
arr=zeros(3,3);
output=uint8(zeros(m,n));
new=zeros(1,9);

for i=2:1:m+1 %starts with 2, and ended with row+1 because after padding, original image inside pad. 
     for j=2:1:n+1 %starts with 2, and ended with column+1 because after padding, original image inside pad. 
         
         arr(1,1)= pad(i-1,j-1);
         arr(2,1)= pad(i,j-1);
         arr(3,1)= pad(i+1,j-1);
         arr(1,2)= pad(i-1,j);
         arr(2,2)= pad(i,j);
         arr(3,2)= pad(i+1,j);
         arr(1,3)= pad(i-1,j+1);
         arr(2,3)= pad(i,j+1);
         arr(3,3)= pad(i+1,j+1);

         x=sort(arr(1,:)); %sort for first row.
         y=sort(arr(2,:)); %sort for second row.
         z=sort(arr(3,:)); %sort for third row.
         new=[x,y,z]; %combine the sorted rows(x,y and z).
         s=sort(new); %sort the new combining array.

         MIN=min(s); %find the min value.
         MAX=max(s); %find the max value.
         output(i-1,j-1)=(MIN+MAX)/2; %find the midpoint filter.
         
     end %used midpoint filter,this filter is best form of uniform filter.
end

figure(4);edge(I);title('edge of noisy4');
figure(5);edge(output);title('edge of recovered4');
figure(6);imshow(output);title('recovered4');

imwrite(output,'recovered4.png');


%% For image-5

clear all;clc;

%read image
I=imread('noisy5.png');
figure(1);imshow(I);title('noisy5.png');

%take the histogram of image
hst=imcrop(I);
figure(2);imshow(hst);
hist=imhist(hst);
figure(3);bar(hist);

%size of image
[m n]=size(I);

a=2;b=10;%we assume the values of a and b.
z=randi(255,m,n);%take the value of z as random.

%form the rayleigh
if (z >= a); pz=(2/b)*(z-a) *exp(-((z-a).^2)/b);

else
    pz=0;

end

ng=pz*sum(hist)/sum(pz);
x=[0:255];hold on;
plot(x,ng,'r-','linewidth',2);

%padding.
pad=uint8(zeros(m+2,n+2));
pad(2:m+1,2:n+1)=I; %starts with 2, because after enlarge with pad our image will start pixel of (2,2).

%create matrices.
arr=zeros(3,3);
output=uint8(zeros(m,n));
new=zeros(1,9);


for i=2:1:m+1 %starts with 2, and ended with row+1 because after padding, original image inside pad. 
     for j=2:1:n+1 %starts with 2, and ended with column+1 because after padding, original image inside pad. 
         
         summ=0;
         arr(1,1)= pad(i-1,j-1);
         arr(2,1)= pad(i,j-1);
         arr(3,1)= pad(i+1,j-1);
         arr(1,2)= pad(i-1,j);
         arr(2,2)= pad(i,j);
         arr(3,2)= pad(i+1,j);
         arr(1,3)= pad(i-1,j+1);
         arr(2,3)= pad(i,j+1);
         arr(3,3)= pad(i+1,j+1);
        
         x=sort(arr(1,:)); %sort for first row.
         y=sort(arr(2,:)); %sort for second row.
         z=sort(arr(3,:)); %sort for third row.
         new=[x,y,z]; %combine the sorted rows(x,y and z).
         s=sort(new); %sort the new combining array.
         MIN=min(s);%find the min value.
         output(i-1,j-1)=MIN; 
         
     end %min filter is useful for finding the darkest points in the image.
end %also, this filter reduces the salt noise from the elephant.

figure(4);edge(I);title('edge of noisy5');
figure(5);edge(output);title('edge of recovered5');
figure(6);imshow(output);title('recovered5');

imwrite(output,'recovered5.png');

%% Bandpass and reject filter
function [hbp,hbr] = BandPassAndRejectFilters(type,m,n,d0,N,w)
%computes frequency domain and bandreject filters.
%h = BandPassAndRejectFilters(type,m,n,d0,N,w) creates a bandpass, hbp
%and a bandreject filter, hbr, of the specified type and size (m,n)

%compute the distances D(u,v)
for u=1:m
    for v=1:n
        D(u,v)=((u-(m/2))^2+(v-(n/2))^2)^(1/2);
    end
end

%compute the filter
switch type
    case 'ideal'
        hbr=ones(m,n);
        hbr(D>=(d0-w/2)& D<=(d0+w/2))=0;
        hbp=1-hbr;
    case 'butterworth'
        hbr=1./(1+(((D.*w)./(D.^2-d0^2)).^(2*N)));
        hbp=1-hbr;
    case 'gaussian'
        hbr=1-exp(-((D.^2-d0^2)./(D*w)).^2);
        hbp=1-hbr;
end
end    

%% Notchpass and reject filter

function [hnp,hnr] = NotchPassAndRejectFilters(type,m,n,d0,uk,vk,N)%computes frequency domain notchpass and notchreject filters.

%h = NotchPassAndRejectFilters(type,m,n,d0,uk,vk,n) creates a notchpass,
%hnp, and a notchreject filter, hnr, of the specified type and size(m,n).

%valid values for type, d0, and n are:

%define D(u,v)
for u=1:m
    for v=1:n
        dkp(u,v)=((u-(m/2)-uk)^2 + (v-(n/2)-vk)^2)^(1/2);
        dkn(u,v)=((u-(m/2)+uk)^2 + (v-(n/2)+vk)^2)^(1/2);
    end
end

%compute the filter
switch type
    case 'ideal' %'ideal' Ideal filter with cutoff frequency d0 and centers(uk,vk).
        h=zeros(m,n);
        h(dkp<=d0)=1;hp=1-h;
        h=zeros(m,n);
        h(dkn<=d0)=1;hn=1-h;
        hnr=hp.*hn;hnp=1-hnr;
    case 'butterworth' %'btw'  Butterworth filter of order n, cutoff d0 and centers(uk,vk).
        h=1./(1+(dkp./d0).^(2*N));hp=1-h;
        h=1./(1+(dkn./d0).^(2*N));hn=1-h;
        hnr=hp.*hn;hnp=1-hnr;
    case 'gaussian' %'gaussian'  Gaussan filter with cutoff frequency d0 and creates centers(uk,vk).
        h=exp(-(dkp.^2)./(2*(d0^2)));hp=1-h;
        h=exp(-(dkn.^2)./(2*(d0^2)));hn=1-h;
        hnr=hp.*hn;hnp=1-hnr;
end
end

