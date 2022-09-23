function [ ret ] = jacobi ()
clc
clear
[A,b,n]=matriz1;
for i = 1 : n
	x0(i)= input('Valor inicial: ');
end
tolerance= input ('Tolerancia: ');
iterations = input('Iteraciones:');
counter = 0;
E = tolerance + 1;
while E > tolerance && counter < iterations
	for i = 1:n
		acum1=0;
		for p = 1 : i - 1
			acum1 = acum1 + A(i,p)*x0(p);
		end
		acum2=0;
		for q = i + 1: n
			acum2 = acum2 + A(i,q)*x0(q);
		end
		x1(i)= (b(i)-acum1 - acum2)/A(i,i);
	end
	E = max(abs(x1-x0));
	x0 = x1;
	counter = counter + 1;
	tabla(counter,1)= counter;
	for i = 1:n
		tabla(counter,i+1)= x0(i);
	end
	tabla(counter,n+2)= E;
end
if E < tolerance
		disp(x0)
		fprintf(' es solución con un error de %g \n',E)
	else
		fprintf('El programa falla con esas iteraciones \n ')
end
disp(tabla)
endfunction
