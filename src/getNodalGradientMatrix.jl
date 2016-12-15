export getNodalGradientMatrix

function getNodalGradientMatrix(mesh::TetraMesh)
# G = getGradientMatrix(mesh)

if isempty(mesh.Grad)

edges = mesh.edges
# number of edges
ne = size(edges, 1)
nn = maximum(edges[:])

# assemble gradient as sparse matrix
i =  [1:ne
      1:ne]

j =  [edges[:, 1]
      edges[:, 2]]

g =   [ ones(ne)
       -ones(ne)]

L = getEdgeLength(mesh)

for k = 1:2*ne
   g[k] *= 1/L[i[k]]    #gradient scaled using edge lengths
end

mesh.Grad = sparse(i, j, g, ne, nn)

end

return mesh.Grad

end
