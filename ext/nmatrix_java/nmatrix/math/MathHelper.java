import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.special.Erf;

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

}