function [ ret ] = gauss_seidel ()
clc
clear
[A,b,n]=matriz1;
for i = 1 : n
	x(i)= input('Valor inicial: ');
end
tolerance= input ('Tolerancia: ');
iterations = input('Iteraciones:');
cont = 0;
E = tolerance + 1;
while E > tolerance && cont < iterations
	aux=x;
	for i = 1:n
		acum1=0;
		for p = 1 : i - 1
			acum1 = acum1 + A(i,p)*x(p);
		end
		acum2=0;
		for q = i + 1: n
			acum2 = acum2 + A(i,q)*x(q);
		end
		x(i)= (b(i)-acum1 - acum2)/A(i,i);
	end
	E = max(abs(x-aux));
	cont = cont + 1;
	tabla(cont,1)= cont;
	for i = 1:n
		tabla(cont,i+1)= x(i);
	end
	tabla(cont,n+2)= E;
end
if E < tolerance
		disp(x)
		fprintf(' es solución con error de %g \n',E)
	else
		fprintf('The program failures in these iterations \n ')
end
disp(tabla)
endfunction
