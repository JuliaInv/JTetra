export getNodalMassMatrix

function getNodalMassMatrix(mesh::TetraMesh, sigma::Array{Float64,1})

	nc = mesh.nc
  nn = mesh.nn

	M = zeros(4,4,nc)

	for i=1:nc
		Ti = mesh.Tetras[i,:]
	  	Xi = mesh.Points[Ti,:]

		M[:,:,i] = getLocalNodalMassMatrix(sigma[i],Xi)

	end
	nodalInd = mesh.Tetras

	ii = zeros(Int64,16*nc)
	jj = zeros(Int64,16*nc)
	# mm = zeros(Float64,16*nc)

	cnt = 1
	for i=1:nc
		for i2=1:4
			for i1=1:4
				ii[cnt] = nodalInd[i,i1]
				jj[cnt] = nodalInd[i,i2]
				# mm[cnt] = M[i1,i2,i]
				cnt    += 1
			end
		end
	end
  mm = vec(M)

	return sparse(ii,jj,mm,nn,nn)
  
end

# Nodal mass matrix is independent of geometry. Store mass matrix for unit coefficient and unit volume.
const M = (1/20) * [
  2 1 1 1
  1 2 1 1
  1 1 2 1
  1 1 1 2]

function getLocalNodalMassMatrix(sigma::Float64,X::Array{Float64})
	v =  getCellVolume(X[1,:], X[2,:], X[3,:], X[4,:])
	return (v * sigma) * M
end
