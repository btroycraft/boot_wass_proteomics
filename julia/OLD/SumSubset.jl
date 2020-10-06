module SumSubset

export sum_subset_close

sum_subset_close = function(mat)

  # For the input matrix, output a function that sums across a given set of indices. Output is a closure, so it contains the required information about the matrix X, but does not copy the actual matrix array.

  return function(subset)

    length_subset = length(subset)
    out = 0.

    for ind1 in 1:length_subset
      for ind2 in (ind1+1):length_subset
        out += mat[subset[ind1], subset[ind2]]
      end
    end

    out /= binomial(length_subset, 2)

    return out
  end
end

end
