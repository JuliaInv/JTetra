export getEdgeLength, getFaceArea, getCellVolume

function getEdgeLength(mesh::TetraMesh)
# get edge length

if isempty(mesh.L)

  ne = mesh.ne
  L = zeros(ne)

  for i=1:ne
    i1 = mesh.edges[i,1]
    i2 = mesh.edges[i,2]

    L[i] = getEdgeLength(mesh.Points[i1,:],  mesh.Points[i2,:])
  end

  mesh.L = L

end

return mesh.L
end


function getFaceArea(mesh::TetraMesh)
# get face size

if isempty(mesh.F)
  nf = mesh.nf
  F = zeros(nf)
  for i=1:nf
    i1 = mesh.faces[i,1]
    i2 = mesh.faces[i,2]
    i3 = mesh.faces[i,3]

    F[i] = getFaceArea(mesh.Points[i1,:], mesh.Points[i2,:], mesh.Points[i3,:])
  end

mesh.F = F
end

return mesh.F
end


function getCellVolume(mesh::TetraMesh)
# get volume

if isempty(mesh.V)

  nc = mesh.nc
  T  = mesh.Tetras
  V = zeros(nc)

  for i=1:nc
    i1 = T[i,1]
    i2 = T[i,2]
    i3 = T[i,3]
    i4 = T[i,4]

    V[i] = getCellVolume(mesh.Points[i1,:], mesh.Points[i2,:], mesh.Points[i3,:], mesh.Points[i4,:])
  end

  mesh.V = V

end

return mesh.V
end


function getEdgeLength(P1::Array{Float64,1}, P2::Array{Float64,1})
	return sqrt((P1[1] - P2[1])^2 + (P1[2] - P2[2])^2 + (P1[3] - P2[3])^2)
end

function getFaceArea(P1::Array{Float64,1}, P2::Array{Float64,1}, P3::Array{Float64,1})
	return norm(cross(P3 - P1, P3 - P2))/2
end

function getCellVolume(P1::Array{Float64,1}, P2::Array{Float64,1}, P3::Array{Float64,1}, P4::Array{Float64,1})
	return 1/6.0*abs(dot(P4 - P3, cross(P4 - P1, P4 - P2)))
end
