function A=MatAlea(N0,Nmax,pct)
%function that allows to sort a matrix of index between N0 and Nmax,
%according to a wished percentage pct of integers in [N0, Nmax]

l=floor(pct*(Nmax-N0)/100);
if l==0
    A=[];
else
    A=(N0+1)+floor(((Nmax-1)-(N0+1)).*rand(1,l));

    %Sort the matrix
    A=sort(A);

    %Test to avoid two (or more) identical index in the matrix A
    for i=2:length(A)
        if A(1,i)==A(1,i-1)
            A(1,i)=A(1,i)+1;
        end
    end
end