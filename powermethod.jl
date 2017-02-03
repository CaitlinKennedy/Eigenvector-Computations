function generateSparceMatrix(rows, columns, density)
  sparcematrix = sprand(rows, columns, density)
  fullmatrix = full(sparcematrix)
    @show fullmatrix #sparcematrix
  return fullmatrix
end

function powerMethod(fullMatrix, iteration)
    x0 = ones(size(fullMatrix, 1) , 1)
    power = fullMatrix * x0
    for i=1:(iteration - 1)
        #@show power #sparse(power)
        power = fullMatrix * power
        power = power / norm(power)
    end
    return power
  #   if length(find(power)) < 1
  #       print()
  #   else
  #       findnz(abs(power))
  #       scale = minimum(findnz(abs(power))[3])
  #     if findfirst(power, scale) == 0
  #         scale = scale * -1
  #     end
  #     #doesn't work if there is a zero value
  #     dominantEigenvector = power / scale
  #     return dominantEigenvector
  #  end
end

function normalizedEig(eigvec)
    eigvec
  for i = 1:size(eigvec, 2)
      #replace with a non zero value
      domEig = sparsevec(eigvec[1:end, i])
      findnz(domEig)
      min = minimum(findnz(abs(domEig))[2])
      if findfirst(domEig, min) == 0
          min = min * -1
      end
      eigvec[1:end, i] = eigvec[1:end, i] / min

  end
    return sparse(eigvec)
end

function inefficientPower(fullMat)
    x =(fullMat)^10 * rand(100,1)
    x = x / norm(x)
    return x/ minimum(findnz(abs(x))[3])
end

fullMat = generateSparceMatrix(100, 100, .01)
fullMat = 1000.*fullMat
#fullMat = [2 -12; 1 -5]
 dominantEigenvector = powerMethod(fullMat, 50 ) #sparse(powerMethod(fullMat, 5 ))
 @show sparse(dominantEigenvector)
#@show inefficientPower(fullMat)
println("\n\n")
allEigenvectors = eigvecs(fullMat)
@show sparse(allEigenvectors)
#@show normalizedEig(eigvecs(fullMat))
