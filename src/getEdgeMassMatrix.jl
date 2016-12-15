export getEdgeMassMatrix

getEdgeMassMatrix(mesh::TetraMesh, sigma::Array{Float64,1}) =
  getEdgeMassMatrix(mesh, [repmat(sigma,1,3) zeros(length(sigma),3)])

function getEdgeMassMatrix(mesh::TetraMesh, Sigma::Array{Float64,2})

	nc = mesh.nc
  ne = mesh.ne
  nn = mesh.nn

	M = zeros(6,6,nc)
   	N = round(Int128,nn)
   	globalEdgeNumbering = arrayToUniqueInt(mesh.edges,N)
   	localEdgeNumbering  = zeros(Int64,nc,6)

	for i=1:nc
		Ti = mesh.Tetras[i,:]
	  	Xi = mesh.Points[Ti,:]

		M[:,:,i] = getLocalEdgeMassMatrix(Sigma[i,:],Xi)

		localEdges = [
      Ti[1] Ti[2]
      Ti[1] Ti[3]
      Ti[1] Ti[4]
      Ti[2] Ti[3]
      Ti[2] Ti[4]
      Ti[3] Ti[4]]
		localEdgeNumbering[i,:] = arrayToUniqueInt(localEdges,N)

	end

	edgeInd = reshape(indexin(localEdgeNumbering,globalEdgeNumbering),nc,6)

	ii = zeros(Int64,36*nc)
	jj = zeros(Int64,36*nc)
  # mm = zeros(Float64,36*nc)

	cnt = 1
	for i=1:nc
  	for i2=1:6
	  	for i1=1:6
				ii[cnt] = edgeInd[i,i1]
				jj[cnt] = edgeInd[i,i2]
        # mm[cnt] = M[i1,i2,i]
				cnt    += 1
			end
		end
	end
  mm = vec(M)

	return sparse(ii,jj,mm,ne,ne)
end

function getLocalEdgeMassMatrix(sigma::Array{Float64},X::Array{Float64})
#              * X4
#             /|\
#            / | \
#           /  |  \
#          /   |   \
#         /    |    \
#     t6 /     |     \t5
#       /      |      \
#      /       |t3     \
#     /        |        \
#    /         |         \
# X3/_________ |t4_______ \ X2
#    .         |         .
#      .       |       .
#        .     |     . t1
#     t2   .   |   .
#            . * .
#             X1

	psi = [ zeros(3) for i = 1:4, j = 1:6 ]
	s = (1 - 1 / sqrt(5)) / 4
	t = 1 - 3 * s

	b1 = X[1,:]
	b2 = X[2,:]
	b3 = X[3,:]
	b4 = X[4,:]
	B  = [b2-b1 b3-b1 b4-b1]
	Binv = inv3X3(B)

	gradLambda1 = Binv' * [-1,-1,-1]
	gradLambda2 = Binv' * [ 1, 0, 0]
	gradLambda3 = Binv' * [ 0, 1, 0]
	gradLambda4 = Binv' * [ 0, 0, 1]

  	L        = norm(b2 - b1)
	psi[1,1] = L * (t * gradLambda2 - s * gradLambda1) #psi_12 = |b2-b1|*psi_12   basis functions
	psi[2,1] = L * (s * gradLambda2 - t * gradLambda1)
	psi[3,1] = L * (s * gradLambda2 - s * gradLambda1)
	psi[4,1] = L * (s * gradLambda2 - s * gradLambda1)

  	L        = norm(b3 - b1)
	psi[1,2] = L * (t * gradLambda3 - s * gradLambda1) #psi_13 = |b3-b1|*psi_13   basis functions
	psi[2,2] = L * (s * gradLambda3 - s * gradLambda1)
	psi[3,2] = L * (s * gradLambda3 - t * gradLambda1)
	psi[4,2] = L * (s * gradLambda3 - s * gradLambda1)

  	L        = norm(b4 - b1)
	psi[1,3] = L * (t * gradLambda4 - s * gradLambda1) #psi_14 = |b4-b1|*psi_14   basis functions
	psi[2,3] = L * (s * gradLambda4 - s * gradLambda1)
	psi[3,3] = L * (s * gradLambda4 - s * gradLambda1)
	psi[4,3] = L * (s * gradLambda4 - t * gradLambda1)

  	L        = norm(b3 - b2)
	psi[1,4] = L * (s * gradLambda3 - s * gradLambda2) #psi_23 = |b3-b2|*psi_23   basis functions
	psi[2,4] = L * (t * gradLambda3 - s * gradLambda2)
	psi[3,4] = L * (s * gradLambda3 - t * gradLambda2)
	psi[4,4] = L * (s * gradLambda3 - s * gradLambda2)

  	L        = norm(b4 - b2)
	psi[1,5] = L * (s * gradLambda4 - s * gradLambda2) #psi_24 = |b4-b2|*psi_24   basis functions
	psi[2,5] = L * (t * gradLambda4 - s * gradLambda2)
	psi[3,5] = L * (s * gradLambda4 - s * gradLambda2)
	psi[4,5] = L * (s * gradLambda4 - t * gradLambda2)

  	L        = norm(b4 - b3)
	psi[1,6] = L * (s * gradLambda4 - s * gradLambda3) #psi_34 = |b4-b3|*psi_34   basis functions
	psi[2,6] = L * (s * gradLambda4 - s * gradLambda3)
	psi[3,6] = L * (t * gradLambda4 - s * gradLambda3)
	psi[4,6] = L * (s * gradLambda4 - t * gradLambda3)

	Sigma = [sigma[1] sigma[4]  sigma[5];
	         sigma[4] sigma[2]  sigma[6];
		 sigma[5] sigma[6]  sigma[3]]

	M = zeros(6,6)
        for q=1:4
	 for k=1:6
          vk = Sigma * psi[q,k]
 	  for j=1:6
            M[j,k] += dot(psi[q,j], vk)
          end
	 end
	end

	v = getCellVolume(X[1,:], X[2,:], X[3,:], X[4,:])
	return (v/4)*M
end

function inv3X3(A)
   Ainv =
	1/det(A) * [
	A[2,2]*A[3,3]-A[2,3]*A[3,2]  A[1,3]*A[3,2]-A[1,2]*A[3,3]  A[1,2]*A[2,3]-A[1,3]*A[2,2];
	A[2,3]*A[3,1]-A[2,1]*A[3,3]  A[1,1]*A[3,3]-A[1,3]*A[3,1]  A[1,3]*A[2,1]-A[1,1]*A[2,3];
	A[2,1]*A[3,2]-A[2,2]*A[3,1]  A[1,2]*A[3,1]-A[1,1]*A[3,2]  A[1,1]*A[2,2]-A[1,2]*A[2,1]
	];
	return Ainv
end
