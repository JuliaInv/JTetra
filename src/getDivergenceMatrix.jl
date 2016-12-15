export getDivergenceMatrix

function getDivergenceMatrix(mesh::TetraMesh)

if isempty(mesh.Div)

  nf = mesh.nf
	nc = mesh.nc
	T  = mesh.Tetras

  ii = zeros(Int64,4*nc)
  jj = zeros(Int64,4*nc)
  kk = zeros(4*nc)

	ii = [1:nc
			  1:nc
			  1:nc
			  1:nc]

	jj = [findRow(mesh.faces,[T[:,1] T[:,2] T[:,3]])
	      findRow(mesh.faces,[T[:,1] T[:,2] T[:,4]])
			  findRow(mesh.faces,[T[:,1] T[:,3] T[:,4]])
			  findRow(mesh.faces,[T[:,2] T[:,3] T[:,4]])]
	kk = [ones(nc)
			 -ones(nc)
			  ones(nc)
			 -ones(nc)]

	#    cnt = 1
	#    for i=1:nc
	#    ii[cnt] = nc
	# 	jj[cnt] = findRow(mesh.faces,[T[i,1] T[i,2] T[i,3]])[1]
	# 	kk[cnt] = 1
	# 	cnt += 1
	#
	#    ii[cnt] = nc
	# 	jj[cnt] = findRow(mesh.faces,[T[i,1] T[i,2] T[i,4]])[1]
	# 	kk[cnt] = -1
	# 	cnt += 1
	#
	#    ii[cnt] = nc
	# 	jj[cnt] = findRow(mesh.faces,[T[i,1] T[i,3] T[i,4]])[1]
	# 	kk[cnt] = 1
	# 	cnt += 1
	#
	#    ii[cnt] = nc
	# 	jj[cnt] = findRow(mesh.faces,[T[i,2] T[i,3] T[i,4]])[1]
	# 	kk[cnt] = -1
	# 	cnt += 1
	# end

  F = getFaceArea(mesh)

  for k = 1:4*nc
    kk[k] *= F[jj[k]]      #divergence scaled using face sizes
  end

	mesh.Div = sparse(ii, jj, kk, nc, nf)

end

return mesh.Div

end
