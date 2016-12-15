function findRow(A::Array,b::Array)
 # finds the first row in A such that A[k,:] = b

 N = round(Int128,maximum(A[:]))
 Av = arrayToUniqueInt(A,N)
 bv = arrayToUniqueInt(b,N)
 return indexin(bv,Av)

end


function arrayToUniqueInt(A::Array{Int64,2},N::Int128)

	k = zeros(Int128,size(A,1))
	for i=1:size(A,1)
		for j=1:size(A,2)
		 	k[i] += A[i,j]*N^(j-1)
		end
	end
	return k
end 
