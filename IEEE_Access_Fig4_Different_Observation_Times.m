% Road Side Unit_Localization_Algorithm

clear
close all
clc

tic

r = [200;0]; % RSU position
R = norm(r).^2;

x = [1;2]; % unknow position to be estimated

fsL = 1; % sampleing frequency (Hz)
fsS = 9;
fsH= 10;

tL = 10; % observation period (s)
tS = 10;
tH = 40;

snr = 30; % SNR in dB

runs = 1000;
 
rmse = zeros(length(tL:tS:tH),length(fsL:fsS:fsH),runs); % LLS
rmse_imp = zeros(length(tL:tS:tH),length(fsL:fsS:fsH),runs); % improved WLLS
Est1 = zeros(length(tL:tS:tH),length(fsL:fsS:fsH),runs);
Est2 = zeros(length(tL:tS:tH),length(fsL:fsS:fsH),runs);
rmseTheory = zeros(length(tL:tS:tH),length(fsL:fsS:fsH),runs); % Theory

itt = 0;
for t = tL:tS:tH
    itt = itt + 1;
    
    % ---------------------------------
    itfs = 0;
    for fs = fsL:fsS:fsH
        itfs = itfs + 1;
        
        num = fs*t;
    
    dt = ones(2,num)*1/fs; % time interval
    % 60 km/h, 16.6667 m/s
    v = [0 6*ones(1,num-1); 0 4*ones(1,num-1)]; % speed m/s
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
        
        crlbN(:,it) = xt(:,it) - r; % for CRLB (42)
        
    end

    % ---------------------------------
    A = zeros(num,3);
    % ---------------------------------
    for it2 = 1:num
        A(it2,:) = [-2*(r-dvcum(:,it2)).' 1]; % Ax = b
    end

%     % ---------------------------------
     snrdB = snr;
    
        for itrun = 1:runs
            sigma2 = data.^2/10^(snrdB/10); % sigma2--square of sigma, here we use: SNR_dB = 10log(d^2/sigma^2)
            noise = sqrt(sigma2).*randn(num,1);
            data_noisy = data + noise;
            b = data_noisy.^2 - R - D + Q; % noisy distance measurements

            % ------- Alg 1: LS ------------
            C = pinv(A'*A);
            xyr = C*A'*b; % LS
            xest = xyr(1:2,1); 

            rmse(itt,itfs,itrun) = norm(xest-x)^2;
%             % ------- Alg 2: Weighted LS ------------
            W = pinv(diag(2*sigma2.^2 + 4*sigma2.*data.^2));
            C2 = pinv(A'*W*A);
            xyr2 = C2*A'*W*b; % LS
            xest2 = xyr2(1:2,1); 
            rmse2(itt,itfs,itrun) = norm(xest2-x)^2;

            % -------     Theoretical Value  ------------
            Cb = pinv(diag(2*sigma2.^2 + 4*sigma2.*data.^2));
            Covb =  pinv(A'*Cb*A);
            rmseTheory(itt,itfs,itrun) = (Covb(1,1)+Covb(2,2));

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
            rmse_imp(itt,itfs,itrun) = norm(x_imp - x)^2;
        end
    end

end

% --------------------  Figure 1 -------------------- 
figure(1)
rmseA = mean(rmse,3); % LS
% rmse2A = mean(rmse2,3); % weighted LS
rmse3A = mean(rmse_imp,3); % improved LS
rmseATheory = mean(rmseTheory,3);

% CRLBA = mean(CRLB,3);
 
ColorSpec1 = {[0 0 1]; % R2014a and Earlier
      [1 0 0];
      [0 0.5 0]; % R2014a and Earlier
      [0 0 1]; % R2014a and Earlier
      [1 0 0];
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
LineStyle = {':o','-v','-.','--x','s','-'};

method = {'  LLS',' WLLS','  (43)'};
num = fsL:fsS:fsH;
strP2 = cell(1,size(length(method)*length(num),1));
for it1 = 1:length(num)
    for it2 = 1:length(method)
        strP2(it2 + (it1-1)*length(method)) =  strcat('f_S = ',mat2str(num(it1)),' Hz, ',method(it2));
    end
    prmse(:,1+length(method)*(it1-1):length(method)*it1) = sqrt([rmseA(:,it1) rmse3A(:,it1) rmseATheory(:,it1)]);
end


pxt = tL:tS:tH;
for it = 1:size(prmse,2)
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
    semilogy(pxt,prmse(:,it),LineStyle{itL},'Color',ColorSpec{itC},'LineWidth',1.5);
 
    hold on
end
 
legend(strP2);

grid

xlabel('Observation time  (s)','fontsize',12)
ylabel('RMSE','fontsize',12)
% axis([tL-10*eps tH+10*eps min(prmse(:))/1.2 max(prmse(:))*1.25])
hold off

 
toc