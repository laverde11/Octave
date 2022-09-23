 function [] = elim_gauss_simple ()
  
  A=matriz1(); 
  b=vectorb();
  n=size(b);
  matriz=[A b];
  matriz
  for i = 1:n-1;
	  for j= i+1:n;
		  mult= (matriz(j,i)/matriz(i,i)); 
			  for k = 1:n+1;
			    matriz(j,k) = matriz(j,k) - mult * matriz(i,k);
			  endfor
        matriz
	  endfor
  endfor
  matriz
endfunction