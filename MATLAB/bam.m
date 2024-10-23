function BAM(infile,outfile)
% This function implements the Blind Acoustic Mask (BAM). In this mask, the
% standard deviation of the noise is estimated in each frame of the noisy
% signal using the DATE estimator. This information is used to calculate a
% parameter dq that represents the proportion of clean signal in each noisy
% frame. This parameter is used to select samples from each frame that are
% preserved to maintain speech intelligibility.
%
% Usage: 
%   BAM(infile,outfile) 
%
% where
%   infile: input .wav filename
%   outfile: output .wav filename

% by: Felipe de Souza Farias and Rosângela Fernandes Coelho
% 10/02/2021
%-------------------------------------------------------------------------



window = 32;
confidence_level = 0.95;
alpha=0.85;
beta=0.15;

[x,Srate] = audioread(infile);



%==================================================================
% Speech Enhancement Parameters
%===================================================================
tam= length(x);
Norms = norm(1,tam);
NoisyData1 = zeros(1,length(x));
Norms = abs(x);
len=floor(window*Srate/1000); % frame size in samples
if rem(len,2)==1, len=len+1; end % force len to be an even number

%===================================================================
%  Noise Detection Parameters
%  (based on the code available at 
%  http://perso.telecom-bretagne.eu/pastor/software)
%===================================================================
dimension=1;
sample_size = len;
h = 1/sqrt((1-confidence_level)*4*sample_size);
kmin = floor((sample_size/2)-h*sample_size);
Psi = sqrt(2)*gamma((1+dimension)/2)/gamma(dimension/2);


j=1;
count=1;

while (j <= tam - len)
    
    observation_norms=Norms(j:j+len-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Added by Zuca (!) SNR for each frame
    
    x_quadro = abs(x(j:j+len-1));
    snr = mean(x_quadro)/std(x_quadro);
    
    if snr == 0
        normalized_threshold = sqrt(dimension);
    else
        normalized_threshold = (1/snr)*acosh(exp(snr^2/2));
    end
    
    % End (Zuca)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    dp_nyq(count) = std(x(j:j+len-1)); %standard deviation of noisy signal in frame
    
    
    %===================================================================
    %  Noise Detection Step
    %===================================================================
    U = sort(reshape(observation_norms,1,numel(observation_norms)));
    possible_estimate = zeros(1,length(U));
    estimate_std = inf;
    for index = 1:length(U)
        possible_estimate(index) = sum(U(1:index))./(index*Psi);
    end
    
    U2 = [U(2:length(U)) inf];
    I = find((U <= normalized_threshold*possible_estimate) & (normalized_threshold*possible_estimate) <= U2);

    if isempty(I) == 0
        J = find(I >= kmin);
        if isempty (J) == 0
            estimate_std  = possible_estimate(I(J(1)));
            ybq(count) = U(I(J(1))); %adicionado_felipe
        else
            estimate_std = possible_estimate(length(possible_estimate));
            ybq(count) = U(length(possible_estimate)); %adicionado_felipe
        end
    end

    
    if isempty(J) && count == 1
        estimate_std1(count) = estimate_std;
    else
        estimate_std1(count)=estimate_std;
    end
    
    
    j=j+len;
    count = count + 1;
    
end

threshold=estimate_std1;
date=mean(estimate_std1);
z=find(estimate_std1>date);
threshold(z)=date;



dq = abs(dp_nyq-threshold); % parameter dq

dq_norm = abs(dp_nyq-threshold)./abs(dp_nyq+threshold); %normalized dq

csi = max(dq_norm,alpha.*threshold); %upper threshold csi



%===================================================================
% Selection of the Noisy Components and
% Speech Signal Reconstruction
%===================================================================

% Defines variable to recover the sign of the amplitude values
sign=x;
z=find(x>=0);
sign(z)=1;
z=find(x<0);
sign(z)=-1;
z=[];

denoi = zeros(tam,1);

j=1;
count=1;
while (j <= tam - len)
    
    w_noisy=Norms(j:j+len-1);
    w_sign=sign(j:j+len-1);
    w_cut(1:len)= threshold(count);
    
    w_denoi=w_noisy-(alpha*w_cut');
    z=find(w_denoi<0);
    w_denoi(z)=beta*w_noisy(z);
    
    w_denoi = w_denoi.*(1-dq_norm(count)).^.5 + w_noisy.*(dq_norm(count)).^.5; %assigns greater importance to frames with lower presence of noise

    w_denoi1=w_sign.*(w_denoi);

    
    denoi(j:j+len-1,1)=w_denoi1;
    
    
    j=j+len;
    count = count + 1;
end

%%%%% last frame
j=j-len;
count = count - 1;

w_noisy=Norms(j:tam);
w_sign=sign(j:tam);
w_denoi = w_noisy-(alpha*threshold(count)'); 
z=find(w_denoi<0);
w_denoi(z)=beta*w_noisy(z); % treats amplitude values lower than 0.
w_denoi1=w_sign.*(w_denoi);

denoi(j:tam,1) = w_denoi1;



audiowrite(outfile,denoi,Srate);



end