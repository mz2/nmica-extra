/*
 * Created on Jul 8, 2004
 */
package biobits.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author thomas
 */
public class MapUtils {
    public static void pushMapping(Map map, Object a, Object b) {
        List l = (List) map.get(a);
        if (l == null) {
            l = new ArrayList();
            map.put(a, l);
        }
        l.add(b);
    }
}
