function generateSparceMatrix(rows, columns, density)
  sparcematrix = sprand(rows, columns, density)
  fullmatrix = full(sparcematrix)
  @show sparcematrix
  return fullmatrix
end

function powerMethod(fullMatrix, iteration)
  x0 = ones(size(fullMatrix, 1) , 1)
  power = fullMatrix * x0
  for i=1:(iteration - 1)
    power = fullMatrix * power
  end
  @show scale = minimum(findnz(power)[3])
  if findfirst(power, scale) == 0
      scale = scale * -1
  end
  #doesn't work if there is a zero value
  dominantEigenvector = power / scale
  return dominantEigenvector
end

function normalizedEig(eigvec)
  for i = 1:size(eigvec, 2)
      #replace with a non zero value
      @show domEig = sparsevec(eigvec[1:end, i])
      @show findnz(domEig)
      @show eigvec[1:end, i] = eigvec[1:end, i] / minimum(findnz(domEig)[3])
  end
  return eigvec[1:end, 1]
end

full = generateSparceMatrix(100, 100, .01)
@show dominantEigenvector = powerMethod(full, 10 )
#@show normalizedEig(eigvecs(full))
