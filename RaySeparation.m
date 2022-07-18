clear all;clc, close all;

nargin=0;
% TEST DATA
    if nargin == 0
        
        fileName = mfilename; fullPath = mfilename('fullpath');
        pathCode = strrep(fullPath, fileName, '');
        if ~isempty(pathCode); cd(pathCode); end
        scrsz = get(0,'ScreenSize'); % get screen size for plotting
        
        % Test 
        path = 'Input';
        
        % in order of "signal strength'
        % ordinary TIFFs, not OME-TIFF   
        file{1} = 'perp_plus_0/p180a.tif'; % green channel
        file{2} = 'perp_plus_45/p135a.tif'; % red channel
        file{3} = 'perp_plus_90/p90a.tif'; % blue channel, only noise
        rgbOrder = [2 1 3];
         
        
        % IMPORT
        maxValue = 2^12 - 1; % input is 12-bit (in 16-bit TIFF though)
                             % taken from .OIB file 
        im = cell(length(file),1);
        for i = 1 : length(file)
            im{i} = imread(fullfile(path, file{i}));
            im{i} = double(im{i}) / maxValue; % scale, 
            
        end
        
        % denoise
        denoisingON = false;
        if denoisingON
            
            matResultsFilename = fullfile(path, 'denoisingResults.mat');
            if exist(matResultsFilename, 'file') == 2

                disp('loading denoising results from disk')
                load(matResultsFilename)
                
            else
            
                try
                    disp('Denoising inputs with BM3D (Anscombe transform), channel: ')
                    % Add Anscombe here transform here, and the inverse
                    for i = 1 : length(im)
                        fprintf('%d ', i)
                        scaleRange = 0.7;  %% ... then set data range in [0.15,0.85], to avoid clipping of extreme values
                        [im_VST, y_sigma, transformLimits] = denoise_anscombeTransform(im{i}, 'forward', scaleRange, []);
                        [NA, denoised_VST] = BM3D(1, im_VST); 
                        [im{i}, ~, ~] = denoise_anscombeTransform(denoised_VST, 'inverse', scaleRange, transformLimits);

                    end
                    fprintf('\n ')
                    save(matResultsFilename, 'im')
                    
                catch err
                    err
                    warning('No BM3D (http://www.cs.tut.fi/~foi/GCF-BM3D/)?')
                    disp('No additional denoising done')
                end
                
            end
        end
        
        % set plot flag
        plotInput = true;
        noOfICs = 3;
        verboseStr = 'on';
                 
    else
        % input arguments 
        verboseStr = 'off';
    end


%% PLOT INPUT
    
        % plot
        if plotInput
            fig = figure('Color', 'w', 'Name', 'Input'); rows = 1; cols = 3;
                set(fig,  'Position', [0.4*scrsz(3) 0.325*scrsz(4) 0.6*scrsz(3) 0.60*scrsz(4)])
            i = 1; sp(i) = subplot(rows,cols,i); imshow(im{i}, []); title(['fa1+ ', num2str(i*0)]);
            i = 2; sp(i) = subplot(rows,cols,i); imshow(im{i}, []); title(['fa1+', num2str(45)]);
            i = 3; sp(i) = subplot(rows,cols,i); imshow(im{i}, []); title(['fa1+ ', num2str(90)]);
        end
  
 %computing temperary parameters fa1 minus fa perpendicular
 fa1_minus_fa_perpendicular=ones(480,640);
 for x=1:480
     for y=1:640
          fa1_minus_fa_perpendicular(x,y)=0.5*atan((im{1}(x,y)+im{3}(x,y)-2*im{2}(x,y))/(im{1}(x,y)-im{3}(x,y))) ;      
     end
 end
 
 %computing I perpentidular and I parallel
 I_perp=ones(480,640);I_para=ones(480,640);
 for x=1:480
     for y=1:640   
         I_perp(x,y)=0.5*(im{1}(x,y)+im{3}(x,y))+0.5*(im{1}(x,y)-im{3}(x,y))/cos(2*fa1_minus_fa_perpendicular(x,y));
         I_para(x,y)=0.5*(im{1}(x,y)+im{3}(x,y))-0.5*(im{1}(x,y)-im{3}(x,y))/cos(2*fa1_minus_fa_perpendicular(x,y));

     end
 end
 
 % define simble varibles
 
 syms Sita;
% Sita=ones(480,640); 
 for i=1:2  %480                %be ware of this: large computation for matlab
     for j=1:640  %640
         Sita(i,j)=sym(['sita',num2str(i),num2str(j)]);
     end
 end
 
 
 %%refractive angle
 syms refra; 
 % based on Snell's law
 k=1.5;  %refraction ratio
 refra=asin(sin(Sita)/k);
 
 %%computting R_s perpendicular and R_s parallel
 
 R_s_perp=(sin(Sita-refra)/sin(Sita+refra)).^2;
 R_s_para=(tan(Sita-refra)/tan(Sita+refra)).^2;
 
 
 %%computting R perpendicular and R parallel
 
 R_perp=2*R_s_perp/(1+R_s_perp);
 R_para=2*R_s_para/(1+R_s_para);
 
 %%%computing L_R and L_B
 syms L_R L_B;
 for x=1:480
     for y=1:640
         L_R(x,y)=2*((1-R_para(x,y))*I_perp(x,y)-(1-R_perp(x,y))*I_para(x,y))/(R_perp(x,y)-R_para(x,y));
         L_B(x,y)=2*(R_para(x,y)*I_perp(x,y)-R_perp(x,y)*I_para(x,y))/(R_para(x,y)-R_perp(x,y));
     end
 end
 
 
 %%%%test the known image for mutual information between L_right and
 %%%%L_background:
 
 L1=fa1_minus_fa_perpendicular;
 L2=fa1_minus_fa_perpendicular;
 
 A=L1;B=L2;
 
 mi(A,A),mi(A,B) 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
     
     
        
        
        
        
        
        
        
        








