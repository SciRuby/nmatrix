//http://stackoverflow.com/questions/709961/determining-if-an-object-is-of-primitive-type

import java.util.*;

public class WrapperType
{
  // How-to-use?
  // public static void main(String[] args)
  // {        Object o = 1;
  //     System.out.println(isWrapperType(String.class));
  //     System.out.println(isWrapperType(o.getClass()));
  // }

  private static final Set<Class<?>> WRAPPER_TYPES = getWrapperTypes();

  public static boolean isWrapperType(Class<?> clazz)
  {
    return WRAPPER_TYPES.contains(clazz);
  }

  private static Set<Class<?>> getWrapperTypes()
  {
    Set<Class<?>> ret = new HashSet<Class<?>>();
    ret.add(Boolean.class);
    ret.add(Character.class);
    ret.add(Byte.class);
    ret.add(Short.class);
    ret.add(Integer.class);
    ret.add(Long.class);
    ret.add(Float.class);
    ret.add(Double.class);
    ret.add(Void.class);
    return ret;
  }
}