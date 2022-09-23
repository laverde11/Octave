function [ ret ] = gauss_seidel_relax ()
clc
clear
[A,b,n]=matriz1;
for i = 1 : n
	x(i)= input('Valor inicial: ');
end
tolerance= input ('Tolerancia: ');
iterations = input('Iteraciones:');
w = input('w:');
counter = 0;
E = tolerance + 1;
while E > tolerance && counter < iterations
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
		x(i)= (((b(i)-acum1 - acum2)/A(i,i))*w) + ((1-w)*x(i));
	end
	E = max(abs(x-aux));
	counter = counter + 1;
	tabla(counter,1)= counter;
	for i = 1:n
		tabla(counter,i+1)= x(i);
	end
	tabla(counter,n+2)= E;
end
if E < tolerance
		disp(x)
		fprintf(' es solucion con error de %g \n',E)
	else
		fprintf('Falla en esas iteraciones \n ')
end
disp(tabla)
endfunction
