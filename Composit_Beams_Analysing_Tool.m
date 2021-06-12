%% Initialize MATLAB
format compact
clear all
close all
clc

%% Define Parameters
E11 = 20010000;
E22 = 1301000;
G12 = 1001000;
V12 = 0.3;
V21 = 0.02;
t = 0.005;


Q11 = E11/(1-V12*V21);
Q22 = E22/(1-V12*V21);
Q12 = V21*E11/(1-V12*V21);
Q66 = G12;
Q16 = 0;
Q26 = 0; 

%% Define Matrixes
matT1 = [45 90 -45 -45 90 45];
matT2 = [0 0 0 0 0 0];
matT3 = [0 45 -45 -45 45 0];
matT4 = [45 -45 90 0 90 0 90 -45 45];
matT5 = [0 90 0 90 0 90 0 90 0];
matT6 = [0 -45 45 0 -45 45 0 -45 45];


matT = matT6; % Select matrix

T = length(matT); 
q11 = zeros(1,T);
q22 = zeros(1,T);
q12 = zeros(1,T);
q16 = zeros(1,T);
q26 = zeros(1,T);
q66 = zeros(1,T);
h = zeros(1,T+1);

for i=1:T+1
    h(i) = (-t*((T/2)-(i-1)));
end

for i=1:T
teta = (matT6(i)*pi/180) ;%Check matrix

q11(i) = Q11*cos(teta).^4+2*(Q12+2*Q66)*sin(teta).^2*cos(teta).^2+Q22*sin(teta).^4;
q22(i) = Q11*sin(teta).^4+2*(Q12+2*Q66)*sin(teta).^2*cos(teta).^2+Q22*cos(teta).^4;
q12(i) = (Q11+Q22-4*Q66)*sin(teta).^2*cos(teta).^2+Q12*(sin(teta).^4+cos(teta).^4);
q66(i) = (Q11+Q22-2*Q12-2*Q66)*sin(teta).^2*cos(teta).^2+Q66*(sin(teta).^4+cos(teta).^4);
q16(i) = (Q11-Q12-2*Q66)*sin(teta)*cos(teta).^3+(Q12-Q22+2*Q66)*sin(teta).^3*cos(teta);
q26(i) = (Q11-Q12-2*Q66)*sin(teta).^3*cos(teta)+(Q12-Q22+2*Q66)*sin(teta)*cos(teta).^3;

end
for i=1:T
    A11 = q11(i)*(h(i+1)-h(i));
end
for i=1:T
    A12 = q12(i)*(h(i+1)-h(i));
end
for i=1:T
    A16 = q16(i)*(h(i+1)-h(i));
end
for i=1:T
    A22 = q22(i)*(h(i+1)-h(i));
end
for i=1:T
    A26 = q26(i)*(h(i+1)-h(i));
end
for i=1:T
    A66 = q66(i)*(h(i+1)-h(i));
end


for i=1:T
    B111 = q11(i)*((h(i+1))^2-(h(i))^2);
end
    B11=0.5*B111;
for i=1:T
    B122 = q12(i)*((h(i+1))^2-(h(i))^2);
end
    B12=0.5*B122;
for i=1:T
    B166 = q16(i)*((h(i+1))^2-(h(i))^2);
end
    B16=0.5*B166;
for i=1:T
    B222 = q22(i)*((h(i+1))^2-(h(i))^2);
end
    B22=0.5*B222;
for i=1:T
    B266 = q26(i)*((h(i+1))^2-(h(i))^2);
end
    B26=0.5*B266;
for i=1:T
    B666 = q66(i)*((h(i+1))^2-(h(i))^2);
end
    B66=0.5*B666;
    
   
for i=1:T
    D111 = q11(i)*((h(i+1))^3-(h(i))^3);
end
    D11=(1/3)*D111;
for i=1:T
    D122 = q12(i)*((h(i+1))^3-(h(i))^3);
end
    D12=(1/3)*D122;
for i=1:T
    D166 = q16(i)*((h(i+1))^3-(h(i))^3);
end
    D16=(1/3)*D166;
for i=1:T
    D222 = q22(i)*((h(i+1))^3-(h(i))^3);
end
    D22=(1/3)*D222;
for i=1:T
    D266 = q26(i)*((h(i+1))^3-(h(i))^3);
end
    D26=(1/3)*D266;
for i=1:T
    D666 = q66(i)*((h(i+1))^3-(h(i))^3);
end
    D66=(1/3)*D666;
    
    A = [A11 A12 A16; A12 A22 A26;A16 A26 A66];
    B = [B11 B12 B16; B12 B22 B26;B16 B26 B66];
    D = [D11 D12 D16; D12 D22 D26;D16 D26 D66];
    
display(matT)
display(A)
display(B)
display(D)