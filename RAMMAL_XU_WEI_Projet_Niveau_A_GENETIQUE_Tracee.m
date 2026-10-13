clear all
clc
clf

NomFichier1 = 'Data_Genetique.txt';
NomFichier2 = 'Tracee_Finale.txt';
M = load ('-ascii' ,NomFichier1) ;
N = load ('-ascii' ,NomFichier2) ;

x=M(:,1);
y=M(:,2);
z=M(:,3);

x2=N(:,1);
y2=N(:,2);

figure(1);
hold on;
grid on;

axis tight;
 plot(x,y);
 plot(x,z);

 plot(x2,y2);


