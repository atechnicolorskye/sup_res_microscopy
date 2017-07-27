function simparm = VB_SPT_simulation_AddPSF_GaussianNoise(data, background, psf,sigma)
% data can be many frames
% simparm.y = pois(data * PSF + background)
% simparm.lambda = data * PSF
% sigma: gaussian noise standard deviation 

[img_size,img_size,N] = size(data); 
simparm.y = zeros(img_size,img_size,N); 

xs = 1:1:img_size; 
xs = xs-0.5;  
[Y,X] = meshgrid(xs,xs); 
R = sqrt(X.^2 + Y.^2);
X = X(:); 
Y = Y(:);
one = 1/(2*pi*(abs(psf).^2)) * exp(-((X-floor(mean(xs))).^2+(Y-floor(mean(xs))).^2)/(2*(psf.^2)));  
one = sum(one(:)); 




for idx = 1:N
    %idx
    
    if mod(idx,100)==1
        fprintf(['building frame ', num2str(idx),'~',  num2str(min(idx+99,N)), '... \n']);
    end
     
   [my_x,my_y] = find(data(:,:,idx)~=0); %high resolution grid
    my_x = my_x - 0.5; 
    my_y = my_y -0.5; 
    true_poses{idx} = [my_y,my_x]; 
    
    if isempty(my_x)
        lambda = zeros(img_size);
    else
        
    one = 1;
    nX = length(my_x); 
    %simparm.nX(idx) = nX; 
    %simparm.idx(idx) = idx; 
    %simparm.molecule_density(idx) = nX/(((img_size-2*bdy/dfactor)*CCD_pitch/1000)^2); 
    
    simparm_i_vec = data(:,:,idx);
    Falw  = simparm_i_vec(find(data(:,:,idx)~=0));  

    %Falw =  scale*ones(1,nX);  
    dfactor = 1; 
    Falx  = [my_x;my_y]/dfactor;  
 
    Xt = X;    
    Yt = Y;

    strall = []; 
    for j = 1:nX
       
      str = [num2str(Falw(j)), '*', num2str(1/one * 1/(2*pi*(abs(psf).^2))), '*', 'exp(-(([', num2str(Xt'),']-x(',num2str(j),')).^2 + ([',num2str(Yt') ,']-x(', num2str(j+nX),')).^2)/', num2str((2*(psf.^2))),')'];
 
    if isempty(strall)
    strall = str; 
    else
    strall = [strall, '+', str];   
    end
 
    end

    func = str2func(['@(x) ', strall]);
    qqq = func(Falx);
    
    lambda = reshape(qqq,img_size,img_size) ; 
    lambda = lambda/sum(lambda(:))*sum(Falw); %normalize lambda
    simparm.lambda(:,:,idx) = lambda; 
    end
    simparm.background = background; 
    %simparm.y(:,:,idx) =  poissrnd(lambda+background); 
    simparm.y(:,:,idx) = normrnd(lambda+background,sigma);
end