close all; clear all
 
%% test on pure data
% simulation
T =   100; 
mu1 = [0,0];
C1 = [3,0;0,4];
Q = [2,0;0,3]; 
%Q = [1,0;0,1]; 
%Q = [7,6;6,7];
A = eye(2); 
%A = [2,3;2,4]; 

%Sigma = [1,0;0,1]; 


%%{
%change Sigma
for i = 1:T
    Sigma(:,:,i) = eye(2)*i ;  %change Sigma to see effect!!!!! smaller sigma, better inference
end
%}

%true
s = zeros(T,2); 
s(1,:) = mvnrnd(mu1,C1);

for i = 2:T
   epsilon = mvnrnd([0,0],Q);
   s(i,:) =  A*s(i-1,:)' + epsilon'; 
end
s_true = s; 

%observations
s_hat = zeros(T,2); 
for i = 1:T
    %s_hat(i,:) = mvnrnd(s(i,:),Sigma);
    s_hat(i,:) = mvnrnd(s(i,:),Sigma(:,:,i));
end

figure; plot(s(:,1),'r'); hold on; plot(s(:,2),'b'); 
plot(s_hat(:,1),'r:'); hold on; plot(s_hat(:,2),'b:'); legend('d-1','d-2','d-1 obs','d-2 obs');

% Calculate H
R = Q\A;  %(inv(Q)*A) 
D1 = -  (  inv(Sigma(:,:,1)) + inv(C1) + A'*inv(Q)*A   ); 
%Dt = -  (  inv(Sigma(:,:,1)) + inv(Q) + A'*inv(Q)*A   ); 
DT =  - (  inv(Sigma(:,:,T)) + inv(Q)  ); 

H = zeros(2*T); 
for i = 1:2:2*T
    if i == 1
        H(i:i+1,i:i+1) = D1;
    elseif i == 2*T-1
        i
        H(i:i+1,i:i+1) =  DT;
    else    
        H(i:i+1,i:i+1) =  -  (  inv(Sigma(:,:,(i+1)/2)) + inv(Q) + A'*inv(Q)*A   ); %Dt
    end
end

for i = 1:T-1
    H((i-1)*2+1:(i-1)*2+2,i*2+1:i*2+2) = R; 
    H(i*2+1:i*2+2,(i-1)*2+1:(i-1)*2+2) = R; 
end


% Calculate G
G1 = mu1*inv(C1) + s_hat(1,:)*inv(Sigma(:,:,1));

%Gt = s_hat(t,:)*inv(Sigma);

G = zeros([2*T,1]); 
for t = 1:T
    if t==1
         G(t*2-1:t*2) =  G1;
    else
        % G(t*2-1:t*2) =  s_hat(t,:)*inv(Sigma);
        G(t*2-1:t*2) =  s_hat(t,:)*inv(Sigma(:,:,t));
    end
end
    
% estimate
s_est_vec = -H\G; 
s_est = reshape(s_est_vec,2,T)'; 

figure; plot(s(:,1),'r','LineWidth',1); hold on; plot(s(:,2),'b','LineWidth',1); 
plot(s_hat(:,1),'r-.','LineWidth',1); hold on; plot(s_hat(:,2),'b-.','LineWidth',1); 
plot(s_est(:,1),'k--','LineWidth',1.5); hold on; plot(s_est(:,2),'k:','LineWidth',1.5); 

legend('d-1 true','d-2 true','d-1 obs','d-2 obs','d-1 estimate','d-2 estimate'); title('Kalman filter')
%[G, s_est_vec]


%isequal(H\G, inv(H)*G)

