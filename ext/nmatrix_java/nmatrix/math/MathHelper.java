import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;

public class MathHelper{

  public static double[] log(double base, double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = FastMath.log(base, arr[i]);
    } 
    return result;
  }

  public static double[] erf(double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = Erf.erf(arr[i]);
    } 
    return result;
  }

  public static double[] erfc(double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = Erf.erfc(arr[i]);
    } 
    return result;
  }

  public static double[] gamma(double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = Gamma.gamma(arr[i]);
    } 
    return result;
  }

  public static double[] round(double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = Math.round(arr[i]);
    } 
    return result;
  }

  public static double[] ldexp(double[] arr1, double[] arr){
    double[] result = new double[arr1.length];
    for(int i = 0; i< arr1.length; i++){
      result[i] = arr1[i] * Math.pow(2, arr[i]);
    } 
    return result;
  }

  public static double[] ldexpScalar(double val, double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = val * Math.pow(2, arr[i]);
    } 
    return result;
  }

  public static double[] ldexpScalar2(double val, double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] = arr[i] * Math.pow(2, val);
    } 
    return result;
  }

  public static double[] hypot(double[] arr1, double[] arr2){
    double[] result = new double[arr1.length];
    for(int i = 0; i< arr1.length; i++){
      result[i] =  Math.sqrt(arr2[i] * arr2[i] + arr1[i] * arr1[i]);
    } 
    return result;
  }

  public static double[] hypotScalar(double val, double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] =  Math.sqrt(arr[i] * arr[i] + val * val);
    } 
    return result;
  }

  public static double[] atan2(double[] arr1, double[] arr2){
    double[] result = new double[arr1.length];
    for(int i = 0; i< arr1.length; i++){
      result[i] =  Math.atan2(arr2[i], arr1[i]);
    } 
    return result;
  }

  public static double[] atan2Scalar(double val, double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] =  Math.atan2(val, arr[i]);
    } 
    return result;
  }

  public static double[] atan2Scalar2(double val, double[] arr){
    double[] result = new double[arr.length];
    for(int i = 0; i< arr.length; i++){
      result[i] =  Math.atan2(arr[i], val);
    } 
    return result;
  }

}