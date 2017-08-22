using MatrixNetworks

function powerMethod(matrix)
    x0 = ones(size(matrix, 1) , 1)
    power = x0 #/ norm(x0)

    for i=1:(100000)
      y = matrix * power
      lam = power' * y
      normr = norm(y - power * lam)
      if normr/abs(lam[1]) < 1.0e-8
        #do something to run again.

        lam[1], normr
        return sparse(y/ norm(y)), dot((matrix * y), y)/dot(y,y)
      end
      power = y / norm(y)

    end
    @printf("Doesnt converge quickly enough\n")
    return [0], norm(power)
  end
 #A + At

mats = matrix_network_datasets()
allEigs = zeros(length(mats), 2)
for name = 1:length(mats)
  matrix = load_matrix_network(mats[name])
  @show mats[name]
  n = size(matrix, 1)
  matrix = matrix + matrix'
  dominantEigenvector, dominantEigenvalue = powerMethod(matrix)
  lam,v = eigs(matrix;nev=3)[1:2]
  if length(dominantEigenvector) > 1
    x = -dominantEigenvector
    x[1] += 1.0
    Q = eye(n) - 2*x*x'/(x'*x)[1]
    B = (Q*matrix*Q')[2:end,2:end]
    eigVecTwo, eigValTwo = powerMethod(B)
    lamB, vB = eigs(B;nev=2)[1:2]
    eigVecTwo = Q*[0.0; eigVecTwo]

    if (1-abs(v[:, 1]'*dominantEigenvector))[1] > .001 || (1-abs(v[:, 2]'*eigVecTwo))[1] > .001
      @printf("\nERROR----ERROR----ERROR----ERROR----ERROR----ERROR   ")
      @show(mats[name])
      @show dominantEigenvector
      @show v[:, 1]
      @printf("\n\n")

    end
    allEigs[name, 1] = dominantEigenvalue
    allEigs[name, 2] = eigValTwo
  end
end
@show allEigs
