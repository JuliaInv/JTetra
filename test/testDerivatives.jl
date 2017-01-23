using JTetra
using Base.Test
include("derivative-test.jl")

println("==== Test derivatives of mass matrices ====")

# read tetrahedral mesh
mesh = getTetraMesh("points.txt","tetras.txt")

# scalar and tensor cell properties for mass matrices
a = exp(randn(mesh.nc))
A = zeros(mesh.nc,6)
for i = 1:mesh.nc
  # generate symmetric positive definite 3x3 matrix
  Ai  = randn(3,3)
  u,V = eig(Ai + Ai')
  B   = V' * diagm(exp(u)) * V
  A[i,1] = B[1,1]
  A[i,2] = B[2,2]
  A[i,3] = B[3,3]
  A[i,4] = B[1,2]
  A[i,5] = B[1,3]
  A[i,6] = B[2,3]
end

# random edge vector
u = rand(mesh.ne)

# derivative test
f(m)    = getEdgeMassMatrix(mesh, m) * u
df(v,m) = getdEdgeMassMatrix(mesh, m, u) * vec(v)

println("Edge mass matrix: isotropic")
checkDerivativeIsotropic,error,order = checkDerivativeMax(f,df,a)
@test checkDerivativeIsotropic

println("Edge mass matrix: anisotropic")
checkDerivativeAnisotropic,error,order = checkDerivativeMax(f,df,A)
@test checkDerivativeAnisotropic

print("Passed!\n")
