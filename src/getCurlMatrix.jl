export getCurlMatrix

function getCurlMatrix(mesh::TetraMesh)
# computes the curl matrix from edges to faces

if isempty(mesh.Curl)

  # incidence edge to face
  faces = mesh.faces
  edges = mesh.edges

  e1 = faces[:, [1, 2]]
  e2 = faces[:, [2, 3]]
  e3 = faces[:, [1, 3]]


  j1 = findRow(edges,e1)
  j2 = findRow(edges,e2)
  j3 = findRow(edges,e3)


  m = size(faces,1)
  n = size(edges,1)
  i = [1:m
       1:m
       1:m]


  j = [j1
       j2
       j3]

  c = [ones(m)
       ones(m)
      -ones(m)]

  L = getEdgeLength(mesh)
  F = getFaceArea(mesh)

  for k = 1:3*m
    c[k] *= (L[j[k]] / F[i[k]])
  end

  mesh.Curl = sparse(i, j, c, m, n)

end

return mesh.Curl

end
