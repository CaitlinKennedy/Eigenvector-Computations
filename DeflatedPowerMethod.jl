function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 #/ norm(x0)

    for i=1:(10000)
      y = matrix * power

      lam = power' * y
      normr = norm(y - power * lam)
      if normr/lam[1] < 1.0e-8
        #do something to run again.
        return sparse(transpose(y/ norm(y))), dot((matrix * y), y)/dot(y,y)
      end
      power = y / norm(y)
    end
    @printf("Doesnt converge quickly enough\n")
    return [0], norm(power)
  end

#function inefficientPower(fullMat)
#    x =(fullMat)^10 * rand(100,1)
#    x = x / norm(x)
#    return x/ minimum(findnz(abs(x))[3])
#end

#@show matrix = sprand(10, 10, .2)
@show matrix = [11 -6 4 -2; 4 1 0 0; -9 9 -6 5; -6 6 -6 7]

dominantEigenvector, dominantEigenvalue = powerMethod(matrix)
@show dominantEigenvector
@show dominantEigenvalue
#@show inefficientPower(fullMat)

println("\n")
if dominantEigenvector != [0]
  x = 1/(dominantEigenvalue*dominantEigenvector[1]) * transpose(matrix[1, :])
  sparse(transpose(x))
  deflatedMatrix = matrix - dominantEigenvalue  * transpose(dominantEigenvector) * x
  eigVecTwo, eigValTwo = powerMethod(deflatedMatrix[2:end, 2:end])
  @show eigVecTwo = transpose(cat(1, [0], transpose(eigVecTwo)))
  eigVecTwo =  (eigValTwo - dominantEigenvalue)*eigVecTwo + dominantEigenvalue * (x * transpose(eigVecTwo))*dominantEigenvector
  @show eigValTwo
  @show eigVecTwo / norm(eigVecTwo)

  println("\n")
  @show eig(full(matrix))[1]
  @show sparse(eig(full(matrix))[2])
  #S@show eigs(matrix)
end
