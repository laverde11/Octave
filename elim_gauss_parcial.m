function [ AB ] = elim_gauss_parcial ()
	A=matriz1;
	b=vectorb;
	AB=[A,b];
  [n]=size(A);
		for k = 1:n-1; 
			mayor = max(abs(AB(k)));
			filamayor=k;
			for s = k+1:n;
				if abs(AB(s,k))>mayor;
					mayor = abs(AB(s,k));
					filamayor=s;
				end
			end
			if ((mayor = 0));
				
			else
				 filamayor ~= k;
					for m = k:n+1;
						aux= AB(filamayor,m);
						AB(filamayor,m)=AB(k,m);
						AB(k,m) = aux;
					endfor
					
			endif
    
			for i = k+1:n;
				multiplicador = AB(i,k)/AB(k,k); 
				for j = k:n+1;
					AB(i,j)= AB(i,j)-multiplicador*AB(k,j);
				endfor 
			endfor 
		endfor
		
endfunction