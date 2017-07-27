function simparm = VB_MPT_Simulation(img_size, scale, psf, T, J, mu1, C1, Q, A, background, noise, noise_gaussian_sigma)

%usesage: 
%  simparm = VB_SPT_simulation(img_size, scale, psf, T, mu1, C1, Q, A, background, 'Poisson', [])
%  simparm = VB_SPT_simulation(img_size, scale, psf, T, mu1, C1, Q, A, background, 'Gaussian', noise_gaussian_sigma)

% simparm.path: the estimated particle xy positions (use the middle of one
% pixel, i.w. 0.5, 1.5 ), transfer to  pixel index ceil(simparm.path)

%{
s1 ~ N(mu1, C1)
s_t ~ N(A*S_{t-1}, Q)

mu1 starts xy coordinates

if noise == 'Poisson' : 
  
% simparm.y = pois(data * PSF + background)
% simparm.lambda = data * PSF

if noise == 'Gaussian' : 
% lambda = data * PSF
% gaussian: mean = lambda+background
% noise_gaussian_sigma: gaussian noise standard deviation 

%}

simparm.psf = psf; 
simparm.T = T; 
simparm.scale = scale; 
simparm.spt_para_name = {'mu1','C1','A','Q'}; 
simparm.spt_para = {mu1, C1, A, Q}; 
simparm.background = background; 
simparm.noise = noise; 
simparm.noise_gaussian_sigma = noise_gaussian_sigma; 
simparm.dfactor = 1; 
%% simulation
simparm.path = zeros(T,2*J);  %size T-by 2J: p1d1 p1d2 p2d1 p2d2 p3d1 p3d2 
simparm.steps = zeros(T-1,2*J);
simparm.path_vec = zeros(2*T*J,1); 

data = zeros(img_size,img_size,T); 

path_min = round(img_size/5); 
path_max = round(img_size*4/5); 


for j=1:J
    
clear s_path    
s = zeros(T,2); 
step = zeros(T,2); 
pre = mu1; 
for i = 1:T
    if i==1
       step(i,:) =  (mvnrnd([0,0],C1)); 
    else
       step(i,:) =  (mvnrnd([0,0],Q)); 
    end
       temp = pre*A + step(i,:); 
       %step(i,:);
       
       %temp(temp<1)=1; 
       %temp(temp>img_size)=img_size; 
       
       temp(temp<path_min)=path_min; 
       temp(temp>path_max)=path_max; 
       
       s(i,:) =  temp ; 
       pre = s(i,:); 
end

 
%scale = 2000;  
%for i = 1:T
%   data(ceil(s(i,1)),ceil(s(i,2)),i) = scale; 
%end
   
%VideoData(data,1)   

s_path = ceil(s); %%ORIGINAL 
%s_path = ceil(s) - 0.5;% mypath 0.5, 1.5,....;  pixel index is ceil(mypath), pixel 1, 2... 

%data = zeros(img_size,img_size,T); 


%scale = 2000;  
for i = 1:T
    data(s_path(i,1),s_path(i,2),i) = data(s_path(i,1),s_path(i,2),i) + scale; 
    %data(ceil(s_path(i,1)),ceil(s_path(i,2)),i) = scale; 
end

for i = 1:T-1
   steps(i,:) = s_path(i+1,:) -  s_path(i,:); 
end

s_path_record = s_path - 0.5; % mypath 0.5, 1.5,....;  pixel index is ceil(mypath), pixel 1, 2... ;
simparm.path(:,2*(j-1)+1:2*j) = s_path_record; 
simparm.steps(:,2*(j-1)+1:2*j) = steps; 
temp = s_path_record'; 
simparm.path_vec(2*T*(j-1)+1:2*T*j) = temp(:); 

%imagesc(data(:,:,1)); drawnow; pause;


end
 



simparm.i = data; 

if strcmp(noise,'Poisson')
sim_obs =  VB_SPT_simulation_AddPSF_PoissonNoise(data, background, psf); 
simparm.y = sim_obs.y; 
simparm.lambda = sim_obs.lambda;
elseif strcmp(noise,'Gaussian')
sim_obs = VB_SPT_simulation_AddPSF_GaussianNoise(data, background, psf, noise_gaussian_sigma); 
simparm.y = sim_obs.y; 
simparm.lambda = sim_obs.lambda;
end