module SumSubset

export sum_subset

sum_subset = function(X)

  # For the input matrix, output a function that sums across a given set of indices. Output is a closure, so it contains the required information about the matrix X, but does not copy the actual matrix array.

  return function(subset)

    length_subset = length(subset)
    total = 0.

    for ind1 in 1:length_subset
      for ind2 in (ind1+1):length_subset
        total += X[subset[ind1], subset[ind2]]
      end
    end

    return total
  end
end

end
