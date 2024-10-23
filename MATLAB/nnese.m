function nnese(infile,outfile,window,confidence_level,snr,alpha,beta)
% This function implements the Non-stationary Noise Estimation for Speech Enhancement (NNESE), 
% proposed in [1]. In the NNESE, the acoustic noise standard deviation is 
% estimated using an adaptation of the d-Dimensional Trimmed Estimator 
% (DATE) on a frame-by-frame basis. The estimated noise standard deviation
% is then subtracted from each corrupted sample to compose the enhanced 
% speech signal.
%
% Usage: 
%   nnese(infile,outfile,window,confidence_level,snr,alpha,beta) 
%
% where
%   infile: input .wav filename
%   outfile: output .wav filename
%   window: frame length in miliseconds (default: 32 ms)
%   confidence_level: parameter of the DATE estimator (default: 95%)
%   snr (rho): minimum signal to noise ratio for the DATE estimator (default: rho = 4)
%   alpha: oversubtraction factor for the speech signal reconstruction (default: 0.35)
%   beta: flooring factor for negative amplitude values (default: 0.65)
%
%   Reference:
%      [1] R. Tavares and R. Coelho, "Speech Enhancement with Nonstationary
%      Acoustic Noise Detection in Time Domain", IEEE Signal Processing
%      Letters, vol. 23, no. 1, pp. 6-10, January 2016.
%-------------------------------------------------------------------------

if nargin ~= 2 && nargin ~= 7
    fprintf('USAGE: nnese(infile,outfile,window,confidence_level,snr,alpha,beta) \n');
    return;
end   

if nargin == 2
   window = 32;
   confidence_level = 0.95;
   snr = 4;
   alpha=0.35;
   beta=0.65;
end

[x,Srate,nbits]=wavread(infile);

%==================================================================
% Speech Enhancement Parameters
%===================================================================
tam= length(x);
Norms = norm(1,tam);
NoisyData1 = zeros(1,length(x));
Norms = abs(x);
len=floor(window*Srate/1000); % frame size in samples
if rem(len,2)==1, len=len+1; end; % force len to be an even number

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
N = dimension/2; a = snr^2/2;    
if snr == 0
    normalized_threshold = sqrt(dimension);
else
    normalized_threshold = (1/snr)*acosh(exp(snr^2/2));
end;

j=1;
count=1;
while (j <= tam - len)

    observation_norms=Norms(j:j+len-1);
    
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
    
    if isempty(I) ==0
        J = find(I >= kmin);
        if isempty (J) == 0
            estimate_std  = possible_estimate(I(J(1)));
        else
            estimate_std = possible_estimate(length(possible_estimate));
        end
    end
    %===================================================================
    
    % plot(possible_estimate);hold on; plot(U,'r');hold on; plot(U2,'b');
    if isempty(J)
       limitante(count)=limitante(count-1);
       estimate_std1(count) = estimate_std1(count-1);
    else
        Cutoff=U(I(J(1)));
        estimate_std1(count)=estimate_std;
        limitante(count)=Cutoff;
    end
    
    j=j+len;
    count = count + 1;
    
end


threshold=estimate_std1;
date=mean(estimate_std1);
z=find(estimate_std1>date);
threshold(z)=date;

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
    w_denoi(z)=beta*w_noisy(z); % treats amplitude values lower than 0.
    w_denoi1=w_sign.*(w_denoi);
    denoi(j:j+len-1,1)=w_denoi1;
    j=j+len;
    count = count + 1;
end

wavwrite(denoi,Srate,nbits,outfile);

end