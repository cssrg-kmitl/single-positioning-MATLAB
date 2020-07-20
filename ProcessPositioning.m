%% Positioning calculation from RINEX 2.11
% Calculate user position by using the single point method
% Original by Napat Tongkasem
% Version 1.00 
% (15/02/2019) - Create the program
% 
% 1. The program need linux command. Cygwin must be installed
% - install Cygwin-setup-x86_64.exe (64-bit ver.)
% or download: http://cygwin.com/install.html
% 
% 2. Main program is ProcessPositioning.m
% 
% 3. We have laboratory website, you can visit
% - http://iono-gnss.kmitl.ac.th/
% =================================================
% CSSRG Laboratory
% School of Telecommunication Engineering
% Faculty of Engineering
% King Mongkut's Institute of Technology Ladkrabang
% Bangkok, Thailand
% =================================================
% Output
% mode 1 = No atmospheric delay correction
% mode 2 = Ionospheric  delay correction
% mode 3 = Tropospheric delay correction
% mode 4 = Tropospheric+Ionospheric delay correction

close all;clear
warning off
tic
% 1. copy RINEX v 2.11 to /RINEX folder
% 2. define RINEX file's name
%% RINEX file
r_o_name = 'KMIT2000.20o'; % observation file's name
r_n_name = 'KMIT2000.20n'; % navigation file's name (if blank, program will downloaded from IGS)
% r_n_name = ''; % navigation file's name

% Setting#1
% =========== Program's path ==========================
p_path = [pwd '\'];             % Program path
R_path = [p_path 'RINEX\'];     % RINEX path
S_path = [p_path 'Results\'];    % Results path
path(path,[p_path 'function']);

%% 1. Read RINEX (using readrinex .mex file)
% Check file
checkfileRN(r_o_name,R_path);
% Read RINEX  
[obs,nav,doy,Year] = readrinex211(r_o_name,r_n_name,R_path); 
year  = num2str(obs.date(1));
month = num2str(obs.date(2),'%.2d');
date  = num2str(obs.date(3),'%.2d');

%% 2. Estimate user positioning
positioning(obs,nav,doy,S_path);


%% 3. Plot Error
ploterror(year,month,date,obs.station,S_path)
toc
% remove file (reset)
% for m = 1:4; rm = num2str(m);delete([S_path 'mode_' rm '\*_' obs.station '_' year '_' month '_' date '.mat']);end
warning on

