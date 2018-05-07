function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 / norm(x0)

    for i=1:(10000)
      y = matrix * power

      lam = power' * y
      newPower = y / norm(y)
      normr = norm(y - power * lam )
      if normr/lam[1] < 1.0e-8
        @show i
        return newPower
      end
      power = newPower
    end
    @printf("Doesnt converge quickly enough\n")
    return zeros(0)
  end

#function inefficientPower(fullMat)
#    x =(fullMat)^10 * rand(100,1)
#    x = x / norm(x)
#    return x/ minimum(findnz(abs(x))[3])
#end

matrix = sprand(1000, 1000, .005) #bigger thqn 3log(n)x/2
matrix = matrix + matrix'
dominantEigenvector = powerMethod(matrix)
dominantEigenvector
#@show inefficientPower(fullMat)
println("\n\n")
if length(dominantEigenvector) == 0
  lam,v = eigs(matrix;nev=2)[1:2] # returns the largest eigenvalue lambda vector v
  @show lam
  @show abs(lam)
else
  lam,v = eigs(matrix;nev=1)[1:2] # returns the largest eigenvalue lambda vector v
  @show abs(lam-(dominantEigenvector'*matrix*dominantEigenvector)[1])
  @show 1-abs(v'*dominantEigenvector) # should be close to 0
end
