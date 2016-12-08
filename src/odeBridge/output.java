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
public class output extends Structure {
    
    public Pointer vals; // double*

    @Override
    protected List getFieldOrder() {
        return Arrays.asList("vals");
    }
    

}
