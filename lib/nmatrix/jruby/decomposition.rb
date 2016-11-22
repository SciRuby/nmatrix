class NMatrix

  # discussion in https://github.com/SciRuby/nmatrix/issues/374

  def matrix_solve rhs
    if rhs.shape[1] > 1
      nmatrix = NMatrix.new :copy
      nmatrix.shape = rhs.shape
      res = []
      #Solve a matrix and store the vectors in a matrix
      (0...rhs.shape[1]).each do |i|
        res << self.solve(rhs.col(i)).s.toArray.to_a
      end
      #res is in col major format
      result = ArrayGenerator.getArrayColMajorDouble res.to_java :double, rhs.shape[0], rhs.shape[1]
      nmatrix.s = ArrayRealVector.new result

      return nmatrix
    else
      return self.solve rhs
    end
  end

end