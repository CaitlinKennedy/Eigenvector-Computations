function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 #/ norm(x0)

    for i=1:(10000)
      y = matrix * power

      lam = power' * y
      newPower = y / norm(y)
      normr = norm(y - newPower * lam)
      if normr/lam[1] < 1.0e-8
        @show i
        return sparse(newPower)
      end
      power = newPower
    end
    @printf("Doesnt converge quickly enough\n")
    return [0]
  end

#function inefficientPower(fullMat)
#    x =(fullMat)^10 * rand(100,1)
#    x = x / norm(x)
#    return x/ minimum(findnz(abs(x))[3])
#end

matrix = sprand(1000, 1000, .0005)
dominantEigenvector = powerMethod(matrix)
@show dominantEigenvector
#@show inefficientPower(fullMat)
println("\n\n")
if dominantEigenvector != [0]
  allEigenvectors = eigvecs(full(matrix))
  @show sparse(allEigenvectors)
end
