public class ArrayComparator{

  public static boolean equals(double[] arr1, double[] arr2){

    double delta = 1e-3;
    
    for(int i=0; i < arr1.length; i++){
      if(Math.abs(arr1[i] - arr2[i]) > delta){
        return false;
      }
    }
    
    return true;
  
  }
}