export getEdgeIntegralOfPolygonalChain

function getEdgeIntegralOfPolygonalChain(mesh::TetraMesh, polygon::Array{Float64,2}; normalize=false)

X     = mesh.Points
T     = mesh.Tetras
edges = mesh.edges
nn    = mesh.nn
ne    = mesh.ne
nc    = mesh.nc

# create connectivity matrix
E = sparse([edges[:,1]; edges[:,2]], [edges[:,2]; edges[:,1]], [edges[:,1]; edges[:,2]], nn, nn)

# number of polygon segments
np = size(polygon,1) - 1

# collect all vertices that discretize the polygon
v  = Int64[]

# locate mesh vertex closest to first point of polygon
X1 = polygon[1,:]
k1 = 0
dmin = Inf
for ii = 1:nn
  d = norm(X[ii,:] - X1)
  if d < dmin
    k1   = ii
    dmin = d
  end
end
X1 = X[k1,:]
push!(v, k1)

# loop over all polygon segments
for ip = 1:np

  # locate mesh vertex closest to end point of current polygon segment
  X2 = polygon[ip+1,:]
  k2 = 0
  dmin = Inf
  for ii = 1:nn
    d = norm(X[ii,:] - X2)
    if d < dmin
      k2   = ii
      dmin = d
    end
  end
  X2 = X[k2,:]

  # find edges connecting k1 (X1) and k2 (X2) by a straight line
  while k1 != k2

    X1 = X[k1,:]
    K  = nonzeros(E[:,k1])  # all vertices connected to k1 by an edge

    # find vertex closest to straight line between k1 (X1) and k2 (X2)
    dmin = Inf
    ll   = norm(X2 - X1)
    for k in K
      Xk = X[k,:]
      l1 = norm(Xk - X1)
      l2 = norm(Xk - X2)
      d  = (l1 + l2) / ll # d should be equal 1 if k (Xi) is on a straight line
      if d < dmin
        k1   = k
        dmin = d
      end
    end
    push!(v, k1)

  end

  # end point of current segment is start point of next segment
  X1 = X2

end

# edge lengths
L = getEdgeLength(mesh)
F = getFaceArea(mesh)
V = getCellVolume(mesh)

# find edge numbers for edges v[1:2], v[2:3], ...
n = length(v) - 1
e = sort(Int64[v[1:n] v[2:n+1]], 2) # sort because edges[:,1] < edges[:,2]
j = findRow(edges, e)

# assemble sparse source vector
s = spzeros(ne,1)
for ii = 1:n
  k = j[ii]
  if v[ii+1] - v[ii] > 0   # line current flows in direction of edge
    s[k] =  L[k]
  else                     # line current flows opposite to direction of edge
    s[k] = -L[k]
  end
end

# normalize
if normalize

	if all(polygon[1,:] .== polygon[np+1,:])

		# closed polygon: divide by enclosed area
		a  = 0.0
		px = polygon[2:np,1] .- polygon[1,1]
		py = polygon[2:np,2] .- polygon[1,2]
		pz = polygon[2:np,3] .- polygon[1,3]
		for ip = 1:np-2
			cx = py[ip] * pz[ip+1] - pz[ip] * py[ip+1]
			cy = pz[ip] * px[ip+1] - px[ip] * pz[ip+1]
			cz = px[ip] * py[ip+1] - py[ip] * px[ip+1]
			a += sqrt(cx * cx + cy * cy + cz * cz)
		end
		a *= 0.5

	else

		# open polygon: divide by length
		a = 0.0
		for ip = 1:np
			dx = polygon[ip+1,1] - polygon[ip,1]
			dy = polygon[ip+1,2] - polygon[ip,2]
			dz = polygon[ip+1,3] - polygon[ip,3]
			a += sqrt(dx * dx + dy * dy + dz * dz)
		end

	end

	s /= a

end

return s

end
