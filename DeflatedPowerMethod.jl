
function powerMethod(matrix, mymatvec)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 #/ norm(x0)

    for i=1:(10000)
      y = matrix * power
      mymatvec += 1
      lam = power' * y
      normr = norm(y - power * lam)
      if normr/abs(lam[1]) < 1.0e-8
        #do something to run again.
        mymatvec += 1
        @show lam[1], normr
        return sparse(y/ norm(y)), dot((matrix * y), y)/dot(y,y), mymatvec
      end
      power = y / norm(y)
    end
    @printf("Doesnt converge quickly enough\n")
    return [0], norm(power), mymatvec
  end
 #A + At
mymatvec = 0
n = 25
matrix = sprand(n, n, .05)
matrix = matrix + matrix'
dominantEigenvector, dominantEigenvalue, mymatvec = powerMethod(matrix, mymatvec)
lam,v = eigs(matrix;nev=3)[1:2] # returns the largest 2 eigenvalues lambda vectors v
a, b, c, d, matvec, e = eigs(matrix; nev=2)
if length(dominantEigenvector) > 1
  x = -dominantEigenvector
  x[1] += 1.0
  Q = eye(n) - 2*x*x'/(x'*x)[1]
  B = (Q*matrix*Q')[2:end,2:end]
  eigVecTwo, eigValTwo, mymatvec = powerMethod(B, mymatvec)
  lamB, vB = eigs(B;nev=2)[1:2] # returns the largest 2 eigenvalues lambda vectors v
  eigVecTwo = Q*[0.0; eigVecTwo]
  @show eigVecTwo
  @show v[:,2]
  @show abs(lam[1]-(dominantEigenvector'*matrix*dominantEigenvector)[1])
  @show 1-abs(v[:, 1]'*dominantEigenvector) # should be close to 0

  print("\n\n")
  #@show abs(lam[2]-(eigVecTwo'*matrix*eigVecTwo)[1])
  #@show 1-abs(v[:, 2]'*eigVecTwo) # should be close to 0
  print("\n EIGS MATVEC = ", matvec)
  print("\n MY MATVEC = ", mymatvec, "\n")

  #@show (lam)
  #@show (lamB)
end
