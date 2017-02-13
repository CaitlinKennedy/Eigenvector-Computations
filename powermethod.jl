function generateSparceMatrix(rows, columns, density)
  sparcematrix = sprand(rows, columns, density)
  fullmatrix = full(sparcematrix)
  return fullmatrix
end

function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0/ norm(x0)

    for i=1:(10000)
      y = matrix * power

      lam = power' * y
      newPower = y / norm(y)
      normr = norm(y - newPower * lam)
      if normr/lam[1] < 1.0e-8 #1.0e-8
        @show i

        return newPower #dot(newPower, power)
      end
      power = newPower
    end
    @printf("Doesnt converge quickly enough\n")
    return [0]
  end

function normalizedEig(matrix)
    @show values, vectors = eigs(matrix)
    @printf("\nHI\n")
  @show   max = maximum(values)
    index = findfirst(values, max)
    return vectors[index,:]

end

function inefficientPower(fullMat)
    x =(fullMat)^10 * rand(100,1)
    x = x / norm(x)
    return x/ minimum(findnz(abs(x))[3])
end

fullMat = generateSparceMatrix(1000, 1000, .0009)
#fullMat = 1000.*fullMat
#fullMat = [2 -12; 1 -5]
 dominantEigenvector = powerMethod(sparse(fullMat)) #sparse(powerMethod(fullMat, 5 ))
 @show sparse(dominantEigenvector)
#@show inefficientPower(fullMat)
println("\n\n")
if dominantEigenvector != [0]
  allEigenvectors = eigvecs(fullMat)
 @show sparse(allEigenvectors)
#  @show normalizedEig(fullMat)

end
