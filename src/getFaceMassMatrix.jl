export getFaceMassMatrix

getFaceMassMatrix(mesh::TetraMesh, sigma::Array{Float64,1}) =
  getFaceMassMatrix(mesh, [repmat(sigma,1,3) zeros(length(sigma),3)])

function getFaceMassMatrix(mesh::TetraMesh, Sigma::Array{Float64,2})

	nc = mesh.nc
  nf = mesh.nf
  nn = mesh.nn

	M = zeros(4,4,nc)
   	N = round(Int128,nn)
   	globalFaceNumbering = arrayToUniqueInt(mesh.faces,N)
   	localFaceNumbering  = zeros(Int64,nc,4)

	for i=1:nc
		Ti = vec(mesh.Tetras[i,:])
	  	Xi = mesh.Points[Ti,:]

		M[:,:,i] = getLocalFaceMassMatrix(Sigma[i,:],Xi)

		localFaces = [
      Ti[1] Ti[2] Ti[3]
      Ti[1] Ti[2] Ti[4]
      Ti[1] Ti[3] Ti[4]
      Ti[2] Ti[3] Ti[4]]
		localFaceNumbering[i,:] = arrayToUniqueInt(localFaces,N)

	end
	faceInd = reshape(indexin(localFaceNumbering,globalFaceNumbering),nc,4)

	ii = zeros(Int64,16*nc)
	jj = zeros(Int64,16*nc)
  # mm = zeros(Float64,16*nc)

	cnt = 1
	for i=1:nc
		for i2=1:4
			for i1=1:4
				ii[cnt] = faceInd[i,i1]
				jj[cnt] = faceInd[i,i2]
        # mm[cnt] = M[i1,i2,i]
				cnt    += 1
			end
		end
	end
  mm = vec(M)

	return sparse(ii,jj,mm,nf,nf)
end




function getLocalFaceMassMatrix(sigma::Array{Float64},X::Array{Float64})
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

	phi = [ zeros(3) for i = 1:4, j = 1:4 ]
	s = (1 - 1 / sqrt(5)) / 4
	t = 1 - 3 * s

	b1  = X[1, :]
	b2  = X[2, :]
	b3  = X[3, :]
	b4  = X[4, :]
  b21 = b2 - b1
  b31 = b3 - b1
  b41 = b4 - b1
	B  = [b21 b31 b41]

  A = getFaceArea(X[2, :], X[3, :], X[4, :])		  # face 234 - opposite vertex b1
	phi[1,4] =  A * B * [s;s;s]		            # phi_4
	phi[2,4] =  A * B * [t;s;s]
	phi[3,4] =  A * B * [s;t;s]
	phi[4,4] =  A * B * [s;s;t]

  A = getFaceArea(X[1, :], X[3, :], X[4, :])	# face 134 - opposite vertex b2
	phi[1,3] = -A * (B * [s;s;s] - b21)	# phi_3
	phi[2,3] = -A * (B * [t;s;s] - b21)
	phi[3,3] = -A * (B * [s;t;s] - b21)
	phi[4,3] = -A * (B * [s;s;t] - b21)

  A = getFaceArea(X[1, :], X[2, :], X[4, :])	# face 124 - opposite vertex b3
	phi[1,2] =  A * (B * [s;s;s] - b31)	# phi_2
	phi[2,2] =  A * (B * [t;s;s] - b31)
	phi[3,2] =  A * (B * [s;t;s] - b31)
	phi[4,2] =  A * (B * [s;s;t] - b31)

  A = getFaceArea(X[1, :], X[2, :], X[3, :])		  # face 123 - opposite vertex b4
	phi[1,1] = -A * (B * [s;s;s] - b41)	# phi_1
	phi[2,1] = -A * (B * [t;s;s] - b41)
	phi[3,1] = -A * (B * [s;t;s] - b41)
	phi[4,1] = -A * (B * [s;s;t] - b41)

	Sigma = [sigma[1] sigma[4]  sigma[5];
	         sigma[4] sigma[2]  sigma[6];
		       sigma[5] sigma[6]  sigma[3]]

   	M = zeros(4,4)
        for q = 1:4
	 for k = 1:4
          vk = Sigma * phi[q,k]
 	  for j = 1:4
            M[j,k] += dot(phi[q,j], vk)
          end
	 end
	end

	v = getCellVolume(X[1,:], X[2,:], X[3,:], X[4,:])
	return (1 / (36 * v)) * M
end
