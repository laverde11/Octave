## Copyright (C) 2015 sgranad1
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## MetodoGaussSeidel

## Author: sgranad1 <sgranad1@PCSAL-48545>
## Created: 2015-04-16

function [ ret ] = MetodoGaussSeidel ()
%Método de Gauss Seidel
%Integrantes:Maria Adelaida Cano, Daniela Jaramillo, Sarita Granada

clc;
clear;
format long; 

disp('METODO DE GAUSS SEIDEL');
disp(' ');
disp('Integrantes:');
disp('Maria Adelaida Cano');
disp('Daniela Jaramillo');
disp('Sarita Granada');
disp(' ');

A=input('Ingresar la matriz de coeficientes  A= ');
b=input('Ingresar la matriz de terminos independientes b= ');
x=input('Ingresar el vector (En columnas) con las aproximaciones iniciales: ');
niter=input('Ingresar el numero de iteraciones: ');
tole=input('Ingresar el valor de la tolerancia: ');

detA=det(A); 		%Se calcula el determinante de la matriz de coeficientes A
if detA==0;			%Verifica que el determinante sea diferente de cero para poder dar solucion a las ecuaciones, si es cero, no tiene unica solucion
disp('El determinanate es cero, el sistema no tiene unica solucion ');
end

n=length(b); 		%Se calcula el numero de elementos del vector b
D=diag(diag(A)); 	%Instrucción para obtener la matriz diagonal D
L=-tril(A,-1); 		%Instrucción para obtener la matriz triangular inferior L
U=-triu(A,1);		%Instrucción para obtener la matriz triangular superior U
fprintf('\n   A=    \n');
disp(' ');
disp(A);
disp(' ');
fprintf('La matriz D calculada es:');
disp(' ');
disp(D);
disp(' ');
fprintf('La matriz triangular inferior L calculada es: \n');
disp(' ');
disp(L);
disp(' ');
fprintf('La matriz triangular superior U calculada es: \n');
disp(' ');
disp(U);
disp(' ');
fprintf('\n Solucion: \n');
fprintf('\n La matriz de transicion de Gauss Seidel es :\n');
Tgs=inv(D-L)*U; 			%Matriz de transición de Gauss Seidel
disp(' ');
disp(Tgs);
respec=max(abs(eig(Tgs))); 	%Calcula el radio espectral para verificar si el método converge
disp(' ');
fprintf(' El radio espectral es igual a: %g',respec);
disp(' ');

if respec>1;
disp(' El radio espectral es mayor que 1 ');
disp(' El metodo no converge ');

	else
	fprintf('\n El vector de iteracion de Gauss Seidel es: \n');
	Cgs=inv(D-L)*b; 					% Vector de iteración Cgs para el método de Jacobi
	disp(' '); 
	disp(Cgs);
	i=0;
	err=tole+1;

	while err>tole & i<niter;		%Garantiza que las soluciones halladas tengan un error pequeño

	xinic=Tgs*x+Cgs;

%disp(xinic)
err=norm(xinic-x); 					%norma 2
%err=max(abs(xinic-x)); %norma 1
%err=norm(xinic-x)/norm(xinic); %norma relativa
x=xinic;

end
fprintf('\n Solucion hallada en %g iteraciones: \n',i);
disp(' ');
for in=1:n;
    fprintf(' X%g=%g\n',in,xinic(in));
end
disp(' '); 
disp(' Con un error de: ');
disp(' ');
fprintf(' %g  \n',err);

end	
end

endfunction


