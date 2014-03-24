clc;

clear;
%{
fid = fopen('./Data/st.txt','rt');
header = textscan(fid,'%s %s %s',1); 
data = textscan(fid,'%s %s %s');
%data = cellfun(@(x) strrep(x,' ','.'),data,'UniformOutput',false);
%clean_data(:,1) = arrayfun(@(x) datenum([data{2}{x}]), 1:length(data{1}) )';
clean_data(:,1:2) = cell2mat(cellfun(@str2num , [data{2:3}],'UniformOutput',false));
fclose(fid);

figure(1);
loglog(clean_data(:,1),clean_data(:,2),'-s');
legend('p=0.6');
ylabel('S(t)');
xlabel('t');

%grid on;

fid = fopen('./Data/fk.txt','rt');
%header2 = textscan(fid,'%s %s %s',1); 
data2 = textscan(fid,'%s %s %s');
%data = cellfun(@(x) strrep(x,' ','.'),data,'UniformOutput',false);
%clean_data(:,1) = arrayfun(@(x) datenum([data{2}{x}]), 1:length(data{1}) )';
clean_data2(:,1) = 1:241;
clean_data2(:,2) = cell2mat(cellfun(@str2num , [data2{1}],'UniformOutput',false));
fclose(fid);

figure(2);
loglog(clean_data2(:,1),clean_data2(:,2),'-s');
legend('k=240');
ylabel('fk');
xlabel('k');
%}
fid = fopen('./Data/areaCCDF0.txt','rt');
header3 = textscan(fid,'#%s',1); 
data3 = textscan(fid,'%s %s');
%data = cellfun(@(x) strrep(x,' ','.'),data,'UniformOutput',false);
%clean_data(:,1) = arrayfun(@(x) datenum([data{2}{x}]), 1:length(data{1}) )';
%clean_data3(:,1) = 1:240;
clean_data3(:,1:2) = cell2mat(cellfun(@str2num , [data3{1:2}],'UniformOutput',false));
fclose(fid);

figure(3);
%x = 1:size(clean_data3(:,1),1);
%X = clean_data3(:,1)';
%Y = clean_data3(:,2)';
loglog(clean_data3(:,1),clean_data3(:,2), '-s');
legend('rho_0=1');
title('human > 0')
ylabel('CCDF');
xlabel('area');


fid = fopen('./Data/densityCCDF.txt','rt');
header4 = textscan(fid,'#%s',1); 
data4 = textscan(fid,'%s %s');
clean_data4(:,1:2) = cell2mat(cellfun(@str2num , [data4{1:2}],'UniformOutput',false));
fclose(fid);
figure(4);
loglog(clean_data4(:,1),clean_data4(:,2), '-s');
legend('rho_0=1');
title('Density')
ylabel('CCDF');
xlabel('density');


fid = fopen('./Data/gridAreaDegreeCCDF.txt','rt');
header5 = textscan(fid,'#%s',1); 
data5 = textscan(fid,'%s %s');
clean_data5(:,1:2) = cell2mat(cellfun(@str2num , [data5{1:2}],'UniformOutput',false));
fclose(fid);
figure(5);
loglog(clean_data5(:,1),clean_data5(:,2), '-s');
legend('rho_0=1');
title('Degree')
ylabel('CCDF');
xlabel('Degree');


fid = fopen('./Data/gridAreaSizeCCDF.txt','rt');
%header6 = textscan(fid,'#%s',1); 
data6 = textscan(fid,'%s %s');
clean_data6(:,1:2) = cell2mat(cellfun(@str2num , [data6{1:2}],'UniformOutput',false));
fclose(fid);
figure(6);
loglog(clean_data6(:,1),clean_data6(:,2), '-s');
legend('rho_0=1');
title('Grid area size');
ylabel('CCDF');
xlabel('Degree');



fid = fopen('./Data/transportProbCCDF.txt','rt');
%header7 = textscan(fid,'#%s',1); 
data7 = textscan(fid,'%s %s');
clean_data7(:,1:2) = cell2mat(cellfun(@str2num , [data7{1:2}],'UniformOutput',false));
fclose(fid);
figure(7);
loglog(clean_data7(:,1),clean_data7(:,2), '-s');
legend('rho_0=1');
title('Prob');
ylabel('CCDF');
xlabel('Prob');
