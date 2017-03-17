function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 #/ norm(x0)

    for i=1:(10000)
      y = matrix * power

      lam = power' * y
      normr = norm(y - power * lam)
      if normr/lam[1] < 1.0e-8
        #do something to run again.
        return sparse(y/ norm(y)), dot((matrix * y), y)/dot(y,y)
      end
      power = y / norm(y)
    end
    @printf("Doesnt converge quickly enough\n")
    return [0], norm(power)
  end

n = 1000
matrix = sprand(n, n, .05)
matrix = matrix + matrix'
dominantEigenvector, dominantEigenvalue = powerMethod(matrix)
lam,v = eigs(matrix;nev=3)[1:2] # returns the largest 2 eigenvalues lambda vectors v

if length(dominantEigenvector) > 1
  x = -dominantEigenvector
  x[1] += 1.0
  Q = eye(n) - 2*x*x'/(x'*x)[1]
  B = (Q*matrix*Q')[2:end,2:end]
  eigVecTwo, eigValTwo = powerMethod(B)
  lamB, vB = eigs(B;nev=2)[1:2] # returns the largest 2 eigenvalues lambda vectors v
  eigVecTwo = Q*[0.0; eigVecTwo]
  @show abs(lam[1]-(dominantEigenvector'*matrix*dominantEigenvector)[1])
  @show 1-abs(v[:, 1]'*dominantEigenvector) # should be close to 0

  print("\n\n")
  @show abs(lam[2]-(eigVecTwo'*matrix*eigVecTwo)[1])
  @show 1-abs(v[:, 2]'*eigVecTwo) # should be close to 0

end
