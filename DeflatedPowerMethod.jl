function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 #/ norm(x0)

    for i=1:(10000)
      y = matrix * power

      lam = power' * y
      newPower = y / norm(y)
      normr = norm(y - newPower * lam)
      if normr/lam[1] < 1.0e-8
        #do something to run again.
        return sparse(newPower), norm(y)
      end
      power = newPower
    end
    @printf("Doesnt converge quickly enough\n")
    return [0], norm(power)
  end

#function inefficientPower(fullMat)
#    x =(fullMat)^10 * rand(100,1)
#    x = x / norm(x)
#    return x/ minimum(findnz(abs(x))[3])
#end

#matrix = sprand(10, 10, .1)
@show matrix = [11 -6  4 -2 ; 4 1 0 0 ; -9 9 -6 5 ; -6 6 -6 7]
dominantEigenvector, dominantEigenvalue = powerMethod(matrix)
@show dominantEigenvector
@show dominantEigenvalue
#@show inefficientPower(fullMat)
@show eigs(matrix)
println("\n\n")
if dominantEigenvector != [0]
  @show matrix[1, :]'
  @show x = 1/(dominantEigenvalue*dominantEigenvector[1]) * matrix[1, :]'
  @show deflatedMatrix = matrix - dominantEigenvalue * dominantEigenvector * x
  @show eigVecTwo, eigValTwo = powerMethod(deflatedMatrix)
  @show (eigValTwo - dominantEigenvalue)*eigVecTwo + dominantEigenvalue * (x * eigVecTwo)*dominantEigenvector
  @show eigValTwo
  allEigenvectors = eigvecs(full(matrix))
  @show sparse(allEigenvectors)
  @show eigs(matrix)[1]
end
