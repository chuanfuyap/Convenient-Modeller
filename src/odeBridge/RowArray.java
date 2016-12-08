/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package odeBridge;

import com.sun.jna.Pointer;
import com.sun.jna.Structure;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author chuanfuyap
 */
public class RowArray extends Structure {
    public Pointer vals;
    public RowArray() {}
		public RowArray(Pointer p) {
			super(p);
		}
    @Override
    protected List getFieldOrder() {
        return Arrays.asList("vals");
    }
    
}
