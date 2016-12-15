using JTetra
using Base.Test

println("==== Test De Rham Complex ====")

mesh = getTetraMesh("points.txt", "tetras.txt")

G = getNodalGradientMatrix(mesh)
C = getCurlMatrix(mesh)
D = getDivergenceMatrix(mesh)

# @printf("||G||     = %.3e\n", norm(G,1))
# @printf("||C||     = %.3e\n", norm(C,1))
# @printf("||D||     = %.3e\n", norm(D,1))
@printf("||C * G|| = %.3e\n", norm(C*G,1))
@printf("||D * C|| = %.3e\n", norm(D*C,1))

@test norm(C*G,1) < 1e-12
@test norm(D*C,1) < 1e-12

print("Passed!\n")
