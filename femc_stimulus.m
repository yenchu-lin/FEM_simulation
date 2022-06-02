function [img,img2,noiseMask] = femc_stimulus(info)
% Create stimulus (letter + noise)
% also see: femc_imge_padding.m
%           femc_Gaussian_noise.m

datafile = cell(info.img.nstim,1);
imgsize = info.img.imsize*info.img.npad;
stim = info.stim;
img = zeros(imgsize, imgsize, info.img.nstim);
img2 = zeros(imgsize, imgsize, info.img.nstim);
if strcmp(info.img.type,'letter')
    datadir = './letter_mat';

    for istim = 1:info.img.nstim

        % Load the input letter.mat files
        datafile{istim} = [datadir filesep stim{istim} '.mat'];
        letter = importdata(datafile{istim},'letter');
        resized_image = imresize(letter,[info.img.imsize,info.img.imsize]);

        % Add padding and adjust contrast of the stimulus
        [I] = femc_imge_padding(resized_image,info.img);

        % Add noise and adjust contrast of the noise
            if strcmp(info.img.noisetype,'Gaussian')
                [img_n_noise,noiseMask] = femc_Gaussian_noise(I,info.img);
            elseif strcmp(info.img.noisetype,'1/F')
                % [img_n_noise,noiseMask] = femc_OneOverF_noise(I,info.img);
                [img_n_noise,noiseMask] = femc_OneOverF_noise2(I,info.img);
            elseif strcmp(info.img.noisetype,'none')
    %             img = I;
    %             img = min(img, 1);
    %             img = max(img, -1);
                img_n_noise = I;
                noiseMask = zeros(imgsize, imgsize);
            end
        img(:,:,istim) = img_n_noise;  
        
    end
    img2 = img;

elseif strcmp(info.img.type,'grating')
    
    for istim = 1:info.img.nstim
        p.freq = info.img.stimsf;       
        if strcmp(stim{istim},'horizontal')
            p.angle(istim) = 90;
        elseif strcmp(stim{istim},'vertical')
            p.angle(istim) = 0;
        end
    end
      x = (0:info.img.imsize-1)*2*pi*(p.freq*info.img.D)/info.img.imsize;  
%     x = (0:imgsize-1)*2*pi*(p.freq*info.img.D*info.img.npad)/imgsize;
    
    % generate sin waves 
    s1 = sin(x); % +1/3*sin(3*x);
    s2 = sin(x); % +1/3*sin(3*x);
%     s1 = cos(x)+1/3*cos(3*x);
%     s2 = sin(x)+1/3*sin(3*x);
    
    % shift of grating
%     shift_ind = round(imgsize/(info.img.stimsf*info.img.D*info.img.npad)*info.img.shift);
%     ss1 = [s1(1,shift_ind+1:end),s1(1,1:shift_ind)];
%     ss2 = [s2(1,shift_ind+1:end),s2(1,1:shift_ind)];

    rs_s1 = (s1-min(s1))/(max(s1)-min(s1)); 
    rs_s2 = (s2-min(s2))/(max(s2)-min(s2));
    adj_rs_s1 = rs_s1 * info.img.cstim;
    adj_rs_s2 = rs_s2 * info.img.cstim;
    adj_s1 = repmat(adj_rs_s1,info.img.imsize,1);
    adj_s2 = repmat(adj_rs_s2,info.img.imsize,1);
    
    % add padding
    n = (imgsize-length(adj_s1))/2; 
    I1 = [zeros(n,imgsize);zeros(length(adj_s1),n),adj_s1,zeros(length(adj_s1),n);zeros(n,imgsize)];
    I2 = [zeros(n,imgsize);zeros(length(adj_s2),n),adj_s2,zeros(length(adj_s2),n);zeros(n,imgsize)];
    
    if strcmp(info.img.noisetype,'Gaussian')
        [img_n_noise1,noiseMask] = femc_Gaussian_noise(I1,info.img);
        [img_n_noise2,~]         = femc_Gaussian_noise(I2,info.img);
    elseif strcmp(info.img.noisetype,'1/F')
        [img_n_noise1,noiseMask] = femc_OneOverF_noise(I1,info.img);
        [img_n_noise2,~]         = femc_OneOverF_noise(I2,info.img);
    elseif strcmp(info.img.noisetype,'none')
        img_n_noise1 = I1;
        img_n_noise2 = I2;
        noiseMask = zeros(imgsize, imgsize);
    end
    J1 = imrotate(img_n_noise1, p.angle(1));
    J2 = imrotate(img_n_noise2, p.angle(2));
    img(:,:,1) = J1;
    img(:,:,2) = J2;
    
    % rotate the image based on the input angles
    img2(:,:,1) = imrotate(I1, p.angle(1));
    img2(:,:,2) = imrotate(I2, p.angle(2));
    
%     for istim = 1:imginfo.nstim       
%         p.sf = imginfo.stimsf(istim);
%         p.width = imginfo.imsize * imginfo.npad; % width of generated image
%         p.height = imginfo.imsize * imginfo.npad; % width of generated image
%         p.px = 0.5; % horizontal center of gabor patch; 0-1
%         p.py = 0.5; % vertical center of gabor patch; 0-1 
%         if strcmp(stim{istim},'left') == 1
%             p.theta = deg2rad(-45); % orientation of the gabor patch; 45 deg right; -45 deg left
%         elseif strcmp(stim{istim},'right') == 1
%             p.theta = deg2rad(45);
%         end
%         p.lambda = p.width/p.sf/(imginfo.npad*imginfo.D); % spatial wavelength
%         p.sigma = 80; % standard deviation of gaussian window
% 
%         [~,~,I] = femc_gabor(p);
%         adjI = I * imginfo.cstim;
%         % datadir = './grating_mat';
%         % datafile{istim} = [datadir filesep stim{istim} '.mat'];
%         % grating = importdata(datafile{istim},'grating');
%         % I = imresize(grating,[imginfo.imsize*imginfo.npad,imginfo.imsize*imginfo.npad]);
%         % [I] = femc_imge_padding(resized_image,imginfo);
%         [img_n_noise] = femc_Gaussian_noise(adjI,imginfo);
%         img(:,:,istim) = img_n_noise;
%     end

end

end