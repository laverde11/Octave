function [ ret ] = Jacobi ()
 
  clc
 
  clear;
  format long; 

  A=input('Ingresar la matriz de coeficientes  A= ');
 
  b=input('Ingresar la matriz de terminos independientes b= ');
 
  x=input('Ingresar el vector (En columnas) con las aproximaciones iniciales: ');

  niter=input('Ingresar el numero de iteraciones: ');
 
  tole=input('Ingresar el valor de la tolerancia: ');

  detA=det(A);
 
    if detA==0;
    disp('El determinanate es cero, el sistema no tiene unica solucion ');
    end
  
  n=length(b);
 
  D=diag(diag(A));
  
  L=-tril(A,-1);
  
  U=-triu(A,1);
  
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
  fprintf('\n La matriz de transicion de Jacobi es :\n');
  
  Tj=inv(D)*(L+U); 
 
  disp(' ');
  disp(Tj);
  
  respec=max(abs(eig(Tj)));
 
  disp(' ');
  fprintf(' El radio espectral es igual a: %g',respec);
  disp(' ');
 
    if respec>1;
    disp(' El radio espectral es mayor que 1 ');
    disp(' El metodo no converge ');
    return
    end
 
  fprintf('\n El vector de iteracion de Jacobi es: \n');

  Cj=inv(D)*b;
 
  disp(' '); 
  disp(Cj);
  i=0;
  err=tole+1;

    while err>tole & i<niter;		
    xinic=Tj*x+Cj;
    err=norm(xinic-x);
    x=xinic;
    i=i+1;
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
endfunction