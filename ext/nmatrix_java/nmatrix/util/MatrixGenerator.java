public class MatrixGenerator
{ 

  // Matrix from Array begin
  public static float[][] getMatrixFloat(float[] array, int row, int col)
  {
    float[][] matrix = new float[row][col];
    for (int index=0, i=0; i < row ; i++){
        for (int j=0; j < col; j++){
            matrix[i][j]= array[index];
            index++;
        }
    }

    return matrix;
     
  }

  public static double[][] getMatrixDouble(double[] array, int row, int col)
  {
    double[][] matrix = new double[row][col];
    for (int index=0, i=0; i < row ; i++){
        for (int j=0; j < col; j++){
            matrix[i][j]= array[index];
            index++;
        }
    }

    return matrix;
     
  }

  public static float[][] getMatrixColMajorFloat(float[] array, int col, int row)
  {
    float[][] matrix = new float[col][row];
    for (int index=0, i=0; i < col ; i++){
        for (int j=0; j < row; j++){
            matrix[i][j]= array[index];
            index++;
        }
    }

    return matrix;
     
  }

  // Matrix from Array end


}