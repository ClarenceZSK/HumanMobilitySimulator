clc;

clear;
%{
fid = fopen('./rho0=1/St.txt','rt');
Stheader1 = textscan(fid,'#%s',1); 
Stdata1 = textscan(fid,'%s %s');
clean_Stdata1(:,1:2) = cell2mat(cellfun(@str2num , [Stdata1{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=2/St.txt','rt');
Stheader2 = textscan(fid,'#%s',1); 
Stdata2 = textscan(fid,'%s %s');
clean_Stdata2(:,1:2) = cell2mat(cellfun(@str2num , [Stdata2{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=3/St.txt','rt');
Stheader3 = textscan(fid,'#%s',1); 
Stdata3 = textscan(fid,'%s %s');
clean_Stdata3(:,1:2) = cell2mat(cellfun(@str2num , [Stdata3{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=4/St.txt','rt');
Stheader4 = textscan(fid,'#%s',1); 
Stdata4 = textscan(fid,'%s %s');
clean_Stdata4(:,1:2) = cell2mat(cellfun(@str2num , [Stdata4{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=5/St.txt','rt');
Stheader5 = textscan(fid,'#%s',1); 
Stdata5 = textscan(fid,'%s %s');
clean_Stdata5(:,1:2) = cell2mat(cellfun(@str2num , [Stdata5{1:2}],'UniformOutput',false));
fclose(fid);

figure(1);
loglog(clean_Stdata1(:,1),clean_Stdata1(:,2), '-s', clean_Stdata2(:,1),clean_Stdata2(:,2), 'b:d', clean_Stdata3(:,1),clean_Stdata3(:,2), 'g-.d', clean_Stdata4(:,1),clean_Stdata4(:,2), 'y-s', clean_Stdata5(:,1),clean_Stdata5(:,2), 'g--p');
legend('rho_0=1', 'rho_0=2', 'rho_0=3', 'rho_0=4', 'rho_0=5');
title('St')
ylabel('S(t)');
xlabel('t');
%}

fid = fopen('./rho0=1/fk.txt','rt');
densitydata1 = textscan(fid,'%s %s');
clean_densitydata1(:,1:2) = cell2mat(cellfun(@str2num , [densitydata1{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=2/fk.txt','rt'); 
densitydata2 = textscan(fid,'%s %s');
clean_densitydata2(:,1:2) = cell2mat(cellfun(@str2num , [densitydata2{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=3/fk.txt','rt');
densitydata3 = textscan(fid,'%s %s');
clean_densitydata3(:,1:2) = cell2mat(cellfun(@str2num , [densitydata3{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=4/fk.txt','rt');
densitydata4 = textscan(fid,'%s %s');
clean_densitydata4(:,1:2) = cell2mat(cellfun(@str2num , [densitydata4{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=5/fk.txt','rt');
densitydata5 = textscan(fid,'%s %s');
clean_densitydata5(:,1:2) = cell2mat(cellfun(@str2num , [densitydata5{1:2}],'UniformOutput',false));
fclose(fid);

figure(2);
loglog(clean_densitydata1(:,1),clean_densitydata1(:,2), '-s', clean_densitydata2(:,1),clean_densitydata2(:,2), 'b:d', clean_densitydata3(:,1),clean_densitydata3(:,2), 'g-.s', clean_densitydata4(:,1),clean_densitydata4(:,2), 'y-d', clean_densitydata5(:,1),clean_densitydata5(:,2), 'g--p');
legend('rho_0=1', 'rho_0=2', 'rho_0=3', 'rho_0=4', 'rho_0=5');
title('FK')
ylabel('f(k)');
xlabel('k');

%{
fid = fopen('./rho0=1/rgt.txt','rt');
degreedata1 = textscan(fid,'%s %s');
clean_degreedata1(:,1:2) = cell2mat(cellfun(@str2num , [degreedata1{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=2/rgt.txt','rt');
degreedata2 = textscan(fid,'%s %s');
clean_degreedata2(:,1:2) = cell2mat(cellfun(@str2num , [degreedata2{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=3/rgt.txt','rt');
degreedata3 = textscan(fid,'%s %s');
clean_degreedata3(:,1:2) = cell2mat(cellfun(@str2num , [degreedata3{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=4/rgt.txt','rt');
degreedata4 = textscan(fid,'%s %s');
clean_degreedata4(:,1:2) = cell2mat(cellfun(@str2num , [degreedata4{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=5/rgt.txt','rt'); 
degreedata5 = textscan(fid,'%s %s');
clean_degreedata5(:,1:2) = cell2mat(cellfun(@str2num , [degreedata5{1:2}],'UniformOutput',false));
fclose(fid);

figure(3);
semilogx(clean_degreedata1(:,1),clean_degreedata1(:,2), '-s', clean_degreedata2(:,1),clean_degreedata2(:,2), 'b:d', clean_degreedata3(:,1),clean_degreedata3(:,2), 'g-.s', clean_degreedata4(:,1),clean_degreedata4(:,2), 'y-d', clean_degreedata5(:,1),clean_degreedata5(:,2), 'g--p');
legend('rho_0=1', 'rho_0=2', 'rho_0=3', 'rho_0=4', 'rho_0=5');
title('RGT')
ylabel('rg(t)');
xlabel('t');


fid = fopen('./rho0=1/Prg.txt','rt');
sizedata1 = textscan(fid,'%s %s');
clean_sizedata1(:,1:2) = cell2mat(cellfun(@str2num , [sizedata1{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=2/Prg.txt','rt');
sizedata2 = textscan(fid,'%s %s');
clean_sizedata2(:,1:2) = cell2mat(cellfun(@str2num , [sizedata2{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=3/Prg.txt','rt');
sizedata3 = textscan(fid,'%s %s');
clean_sizedata3(:,1:2) = cell2mat(cellfun(@str2num , [sizedata3{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=4/Prg.txt','rt');
sizedata4 = textscan(fid,'%s %s');
clean_sizedata4(:,1:2) = cell2mat(cellfun(@str2num , [sizedata4{1:2}],'UniformOutput',false));
fclose(fid);

fid = fopen('./rho0=5/Prg.txt','rt');
sizedata5 = textscan(fid,'%s %s');
clean_sizedata5(:,1:2) = cell2mat(cellfun(@str2num , [sizedata5{1:2}],'UniformOutput',false));
fclose(fid);

figure(4);
loglog(clean_sizedata1(:,1),clean_sizedata1(:,2), '-s', clean_sizedata2(:,1),clean_sizedata2(:,2), 'b:d', clean_sizedata3(:,1),clean_sizedata3(:,2), 'g-.s', clean_sizedata4(:,1),clean_sizedata4(:,2), 'y-d', clean_sizedata5(:,1),clean_sizedata5(:,2), 'g--p');
legend('rho_0=1', 'rho_0=2', 'rho_0=3', 'rho_0=4', 'rho_0=5');
title('Distribution of Rg');
ylabel('P(rg)');
xlabel('rg');
%}



