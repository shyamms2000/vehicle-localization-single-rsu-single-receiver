% Road Side Unit_Localization_Algorithm
clear
close all
clc

tic

r = [200;0]; % RSU position
R = norm(r).^2;

x = [1;2]; % unknow position to be estimated

fs = 10; % sampleing frequency (Hz)

tL = 30; % observation period (s)
tS = 10;
tH = 30;

numL = tL*fs; % number of measurements
numS = tS*fs;
numH = tH*fs;

snrL = 10; % SNR in dB
snrS = 2;
snrH = 30;

runs = 1000;

% ---------------------------------
rmse = zeros(length(numL:numS:numH),length(snrL:snrS:snrH),runs); % LS
rmse2 = zeros(length(numL:numS:numH),length(snrL:snrS:snrH),runs); % WLS
rmse_imp = zeros(length(numL:numS:numH),length(snrL:snrS:snrH),runs); % improved WLS
Est1 = zeros(length(numL:numS:numH),length(snrL:snrS:snrH),runs);
Est2 = zeros(length(numL:numS:numH),length(snrL:snrS:snrH),runs);
rmseTheory = zeros(length(numL:numS:numH),length(snrL:snrS:snrH),runs); % Theory

itnum = 0;
for num = numL:numS:numH
    itnum = itnum + 1; 
    
    dt = ones(2,num)*1/fs; % time interval
    % 60 km/h, 16.6667 m/s
    v = [0 6*ones(1,num-1); 0 4*ones(1,num-1)]; % speed m/s
    dv = v.*dt; % accumulate
    % ---------------------------------
    xt = zeros(2,num);
    dvcum = zeros(2,num);
    data = zeros(num,1);
%     data2 = zeros(num,1);
    D = zeros(num,1);
    Q = zeros(num,1);
    % ---------------------------------

    for it = 1:num
        dvcum(:,it) = sum(dv(:,1:it),2);
        xt(:,it) = x + dvcum(:,it); % real time position
        data(it,1) = norm(xt(:,it) - r); % measured distance m (noiseless)    
        D(it,1) = norm(dvcum(:,it))^2;
        Q(it,1) = 2*dvcum(:,it).'*r;
        
        crlbN(:,it) = xt(:,it) - r; % for CRLB (42)
        
    end

    % ---------------------------------
    A = zeros(num,3);
    % ---------------------------------
    for it2 = 1:num
        A(it2,:) = [-2*(r-dvcum(:,it2)).' 1]; % Ax = b
    end

    % ---------------------------------
    itsnr = 0;
    for snrdB = snrL:snrS:snrH
        itsnr = itsnr + 1;
    
        for itrun = 1:runs
            sigma2 = data.^2/10^(snrdB/10); % sigma2--square of sigma, here we use: SNR_dB = 10log(d^2/sigma^2)
            noise = sqrt(sigma2).*randn(num,1);
            data_noisy = data + noise;
            b = data_noisy.^2 - R - D + Q; % noisy distance measurements

            % ------- Alg 1: LS ------------
            C = pinv(A'*A);
            xyr = C*A'*b; % LS
            xest = xyr(1:2,1); 

            rmse(itnum,itsnr,itrun) = norm(xest-x)^2;
%             % ------- Alg 2: Weighted LS ------------
            W = pinv(diag(2*sigma2.^2 + 4*sigma2.*data.^2));
            C2 = pinv(A'*W*A);
            xyr2 = C2*A'*W*b; % LS
            xest2 = xyr2(1:2,1); 
            rmse2(itnum,itsnr,itrun) = norm(xest2-x)^2;

            % -------     Theoretical Value  ------------
            Cb = pinv(diag(2*sigma2.^2 + 4*sigma2.*data.^2));
            Covb =  pinv(A'*Cb*A);
            rmseTheory(itnum,itsnr,itrun) = (Covb(1,1)+Covb(2,2));
            
            % --------------  CRLB Value ----------------- 
            crlbD = (sigma2.^2 + 2*sigma2.*data.^2).'*2;
            for c1 = 1:2
                for c2 = 1:2
                    crlb2(c1,c2) = sum(4*crlbN(c1,:).*crlbN(c2,:)./crlbD);
                end
            end
            CRLB2 = inv(crlb2);
            CRLB(itnum,itsnr,itrun) = trace(CRLB2);
            
            % -------------------------------------------------------------

% %         % second step
            z = xyr2;
            Z(3) = z(1)^2 + z(2)^2;
            s = sign(z(1:2));
            G = [1 0;0 1;1 1];
            h = [z(1)^2;z(2)^2;z(3)];
            Phi = pinv(diag([2*z(1:2);1])*Covb*diag([2*z(1:2);1]));
            z = pinv(G'*Phi*G)*G'*Phi*h;
            x_imp = real(sign(s).*sqrt(z));
            rmse_imp(itnum,itsnr,itrun) = norm(x_imp - x)^2;
        end
    end

end

% --------------------  Figure 1 -------------------- 
figure(1)
rmseA = mean(rmse,3); % LLS
rmse3A = mean(rmse_imp,3); % improved WLLS
rmseATheory = mean(rmseTheory,3);

CRLBA = mean(CRLB,3);
 
ColorSpec1 = {[0 0 1]; % R2014a and Earlier
      [1 0 0];
      [0 0 0];
      [0 0.5 0]; % R2014a and Earlier
      [0 0 1]; % R2014a and Earlier
      [1 0 0];
      [0 0.5 0];
      [0 0 0]};
%       [0 0.75 0.75];
%       [0.75 0 0.75];
%       [0.75 0.75 0];
%       [0.25 0.25 0.25]};
  
ColorSpec2 = {[0    0.4470    0.7410]; % Starting in R2014b
    [0.8500    0.3250    0.0980];
    [0.9290    0.6940    0.1250];
    [0.4940    0.1840    0.5560];
    [0.4660    0.6740    0.1880];
    [0.3010    0.7450    0.9330];
    [0.6350    0.0780    0.1840]};

ColorSpec = ColorSpec1;
LineStyle = {'-.o','-s',':','x','--*','s','-','x'};

sr = [numL:numS:numH numL:numS:numH numL:numS:numH numL:numS:numH];
method = {'  LLS',' WLLS','  (43)',' CRLB'};
num = numL:numS:numH;
strP2 = cell(1,size(length(method)*length(num),1));
for it1 = 1:length(num)
    for it2 = 1:length(method)
        strP2(it2 + (it1-1)*length(method)) =  strcat('T_o= ',mat2str(num(it1)/fs),'s, ',method(it2));
    end
    prmse(1+length(method)*(it1-1):length(method)*it1,:) = sqrt([rmseA(it1,:);rmse3A(it1,:);rmseATheory(it1,:);CRLBA(it1,:)]);
end


psnr = snrL:snrS:snrH;
for it = 1:size(prmse,1)
    % LineStyle
    if mod(it,length(LineStyle))    
       itL = mod(it,length(LineStyle));
    else
       itL = length(LineStyle);
    end
    % ColorSpec
    if mod(it,length(ColorSpec))    
       itC = mod(it,length(ColorSpec));
    else
       itC = length(ColorSpec);
    end
    semilogy(psnr,prmse(it,:),LineStyle{itL},'Color',ColorSpec{itC},'LineWidth',2);
    hold on
end

legend(strP2);

grid

xlabel('SNR (dB)','fontsize',12)
ylabel('RMSE','fontsize',12)
axis([snrL-10*eps snrH+10*eps min(prmse(:))/1.2 max(prmse(:))*1.25])
hold off

toc