%% SPT 
% single particle

img_size=32; T=100; scale = 2000 %500%8000%4000;  %change scale to change SNR 

mu1 = [round(img_size/2),round(img_size/2)];
C1 = [3,0;0,3];
Q = C1;
A = eye(2); 


psf = 1;
background = 1000; %constant background included into mean (lambda) of the noise 
offset = 0;     %y = poisson(lambda+background) + offset; 
dfactor=1; 


noise =   'Poisson' %
noise_gaussian_sigma = []

 %{
noise =  'Gaussian'
noise_gaussian_sigma =  100  % standard deviation, not variance of gaussian

%}
%% simlations
% true: simparm.path  (true is the index of pixel containning the particle. note that matlab index starts from 1, not 0)
% observations: simparm.y
simparm = VB_SPT_simulation(img_size, scale, psf, T, mu1, C1, Q, A, background, noise, noise_gaussian_sigma);


%write mat file to tiff file, can read in to python 
filename = VB_WriteImage(simparm.y, './', 'obs.tif', T);  

%% see a video
close all
VideoData(simparm.y,1)