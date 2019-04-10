% Road Side Unit_Localization_Algorithm

clear
close all
clc

tic

r = [200;0]; % RSU position
R = norm(r).^2;

x = [1;2]; % unknow position to be estimated

t = 30; % observation period (s)

distanceBias = [4 3 1]; % distance with error, mean
distanceStd = [2 2 1]; % distance with error, standard deviation

fs = 10;

snrL = 20; % SNR in dB
snrS = 2;
snrH = 30;

runs = 1000; % 500runs/30mins

% ---------------------------------
rmse = zeros(length(fs),length(snrL:snrS:snrH),runs);
rmse2 = zeros(length(fs),length(snrL:snrS:snrH),runs);
rmse_imp = zeros(length(fs),length(snrL:snrS:snrH),runs);
Est1 = zeros(length(fs),length(snrL:snrS:snrH),runs);
Est2 = zeros(length(fs),length(snrL:snrS:snrH),runs);


for itnum = 1:length(distanceStd)

    disp(itnum)
    num = t*fs;
    dt = ones(2,num)*1/fs; % time interval
    % 60 km/h, 16.6667 m/s
    v = [0 6*ones(1,num-1); 0 8*ones(1,num-1)]; % speed m/s
    dv = v.*dt; % accumulate
    % ---------------------------------
    xt = zeros(2,num);
    dvcum = zeros(2,num);
    data = zeros(num,1);
    D = zeros(num,1);
    Q = zeros(num,1);
    % ---------------------------------

    for it = 1:num
        dvcum(:,it) = sum(dv(:,1:it),2);
        xt(:,it) = x + dvcum(:,it); % real time position
        data(it,1) = norm(xt(:,it) - r); % measured distance m (noiseless)    
        D(it,1) = norm(dvcum(:,it))^2;
        Q(it,1) = 2*dvcum(:,it).'*r;
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
%             sigma2 = data/10^(snrdB/10); % sigma2--square of sigma, here we use: SNR_dB = 10log(d/sigma^2)
            noise = distanceBias(itnum)*ones(num,1) + distanceStd(itnum)*randn(num,1) + sqrt(sigma2).*randn(num,1);
            data_noisy = data + noise;
            b = data_noisy.^2 - R - D + Q; % noisy distance measurements
            
            % ------- Alg 1: LS ------------
            C = pinv(A'*A);
            xyr = C*A'*b; % LS
            xest = xyr(1:2,1); 
            rmse(itnum,itsnr,itrun) = norm(xest-x)^2;
            
%             % ------- Alg 2: Weighted LS ------------
             W = pinv(diag(2*(distanceStd(itnum)^2+sigma2).^2 + 4*(distanceStd(itnum)^2+sigma2).*data.^2));
             C2 = pinv(A'*W*A);
             xyr2 = C2*A'*W*b; % LS
             xest2 = xyr2(1:2,1); 
             rmse2(itnum,itsnr,itrun) = norm(xest2-x)^2;

%             % -------     Theoretical Value  ------------
            Cb = pinv(diag(2*(distanceStd(itnum)^2+sigma2).^2 + 4*(distanceStd(itnum)^2+sigma2).*data.^2));
            Covb =  pinv(A'*Cb*A);
% %         % second step
             z = xyr2;
             s = sign(z(1:2));
             G = [1 0;0 1;1 1];
             h = [z(1)^2;z(2)^2;z(3)];
             Phi = pinv(diag([2*z(1:2);1])*Covb*diag([2*z(1:2);1]));
             z = pinv(G'*Phi*G)*G'*Phi*h;
             x_imp = real(sign(s).*sqrt(z));
             rmse_imp(itnum,itsnr,itrun) = norm(x_imp - x)^2;
             Est1(itnum,itsnr,itrun) = x_imp(1)-x(1);
             Est2(itnum,itsnr,itrun) = x_imp(2)-x(2);
        end
    end

end

% --------------------  Figure 1 -------------------- 
figure(1)
rmseA = mean(rmse,3); % LLS
rmse3A = mean(rmse_imp,3); % improved WLLS

ColorSpec1 = {[0 0 1]; % R2014a and Earlier
      [1 0 0];[0 0 1]; % R2014a and Earlier
      [1 0 0];[0 0 1]; % R2014a and Earlier
      [1 0 0];[0 0 1]; % R2014a and Earlier
      [1 0 0];
      [0 0 0];
      [0 0 1]; % R2014a and Earlier
      [1 0 0];[0 0.5 0]};
%       [0 0.75 0.75];
%       [0.75 0 0.75];
%       [0.75 0.75 0];
%       [0.25 0.25 0.25]};
  
ColorSpec2 = {[0    0.4470    0.7410]; % Starting in R2014b
    [0.9290    0.6940    0.1250];
    [0.8500    0.3250    0.0980]};

ColorSpec = ColorSpec1;
LineStyle = {':x','-v',':o','-s',':s','-o','-.+',':s','-'};

sr = [distanceStd distanceStd distanceStd distanceStd];
method = {' LLS',' WLLS'};
strP2 = cell(1,size(length(method)*length(distanceStd),1));
for it1 = 1:length(distanceStd)
    for it2 = 1:length(method)
        strP2(it2 + (it1-1)*length(method)) =  strcat('Bias =',' ',num2str(distanceBias(it1)),', Std =',' ',num2str(distanceStd(it1)),', ',method(it2));
    end
    prmse(1+length(method)*(it1-1):length(method)*it1,:) = sqrt([rmseA(it1,:);rmse3A(it1,:)]);
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
    plot(psnr,prmse(it,:),LineStyle{itL},'Color',ColorSpec{itC},'LineWidth',1.5);

    hold on
end

legend(strP2);

grid

xlabel('SNR (dB)','fontsize',12)
ylabel('RMSE','fontsize',12)
xlim([snrL-0.01 snrH])
ylim([1 7])
set(gca,'Xtick',snrL:snrS:snrH)
hold off

toc