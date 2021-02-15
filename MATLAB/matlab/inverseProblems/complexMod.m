clc
clear all

nunits = 1;
freqs = [100:100:1000];
omega = 2*pi*freqs;

c = [1];
d = [1];
tau = [1e-3];

G0 = 20e9;
Ginf = 10E9;
K0 = 200e9;
Kinf = 100e9;


Greal(1:length(freqs))=Ginf;
Gim(1:length(freqs))=0;
Kreal(1:length(freqs))=Kinf;
Kim(1:length(freqs))=0;

for j=1:length(freqs)
    for i=1:nunits
           Greal(j) = Greal(j) + (G0-Ginf) * c(i) * omega(j)^2 * tau(i)^2 / (1 + (omega(j)* tau(i))^2 );
           Gim(j) = Gim(j) + (G0-Ginf) * c(i) * omega(j) * tau(i)/ (1 + (omega(j) * tau(i))^2 );
           Kreal(j) = Kreal(j) + (K0-Kinf) * c(i) * omega(j)^2 * tau(i)^2 / (1 + (omega(j)* tau(i))^2 );
           Kim(j) = Kim(j) + (K0-Kinf) * c(i) * omega(j) * tau(i)/ (1 + (omega(j) * tau(i))^2 ); 
    end
end
Geq = sqrt( Greal.^2 + Gim.^2);
tandelta = Gim ./ Greal;
Keq = sqrt( Kreal.^2 + Kim.^2);
tandeltaK = Kim ./ Kreal;

figure(1)
plot(freqs, Geq);
figure(2)
plot(freqs, tandelta, '--');
figure(3)
plot(freqs, Keq);
figure(4)
plot(freqs, tandeltaK, '--');
           
fid_gre = fopen('Greal.txt','w');
fid_gim = fopen('Gim.txt','w');
fid_kre = fopen('Kreal.txt','w');
fid_kim = fopen('Kim.txt','w');


for p=1:length(freqs)
    fprintf(fid_gre, 'data %g %10.9f  \n', freqs(p), Greal(p));
    fprintf(fid_gim, 'data %g %10.9f \n', freqs(p), Gim(p));
    fprintf(fid_kre, 'data %g %10.9f \n', freqs(p), Kreal(p));
    fprintf(fid_kim, 'data %g %10.9f \n', freqs(p), Kim(p));
end
