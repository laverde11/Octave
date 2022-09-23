function [ ret ] = MetodoGaussianaSimple ()

clc;

clear;
format long;

A=input('Ingrese la matriz de coeficientes  A= ');
b=input('Ingrese la matriz de terminos independientes  b= ');

[n,m]=size(A);
Ab=[A,b];

fprintf('\n La matriz aumentada es Ab= \n');
disp(' ');
disp(Ab);
disp(' ');

	if n==m
   
		for k=1:(n-1)
		   disp(' ');    
     
       fprintf('\n\n ETAPA %g= ',k);
		   disp(' ');            
		  
			 fprintf('Los multiplicadores para esta etapa son: \n');			
           for i=(k+1):n
                   
              m(i,k)= Ab(i,k)/Ab(k,k); 
             
			        fprintf('\n m(%g,%g)= ',i,k)              
			        disp(m(i,k));     
				         for j=k:(n+1)
                  
                   Ab(i,j)= Ab(i,j) - m(i,k)*Ab(k,j); 
				         end
			    end
        
			  fprintf('\n La matriz resultante al final de esta estapa es  Ab= \n');
			  disp(' ');
			  disp(Ab);
			  disp(' ');
		    end	
		 
		  for i=n:-1:1
        
         suma=0;
        
			     for p=(i+1):n
          
             suma = suma + Ab(i,p)*X(p);
			     end
      
         X(i)=(Ab(i,n+1)-suma)/Ab(i,i);
		  end
	
		disp(' ');
		fprintf('\n\n La matriz final es: Ab= \n');
		disp(' ');
		disp(Ab);
		disp(' ');
	
		fprintf('\n\n La solucion (X1,X2,...,Xn) es: \n');
		disp(' ');
      
		  for i=1:n
		  Xi=X(1,i);
		  fprintf('\n X%g=',i)
		  disp(Xi);
		  end
   
	else 
     fprintf('ERROR: LA MATRIZ NO ES CUADRADA');
	   disp(' '); 
	end
	
endfunction