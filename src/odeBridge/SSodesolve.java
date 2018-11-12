/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package odeBridge;

import com.sun.jna.Library;
import com.sun.jna.Memory;
import com.sun.jna.Native;
import com.sun.jna.Pointer;

/**
 *
 * @author chuanfuyap
 */
public class SSodesolve {
    private String link;
    private int parametercount;
    private String[] varyingmetabolites;
    private String[] parameterNames;
    private String[] reactionID;
    private double[] outputarray;
    
    
    public SSodesolve(String sbmllink, String[] variables, int count, String[] names, String[] rID){
        link=sbmllink;
        parametercount=count;
        varyingmetabolites=variables;
        parameterNames=names;
        reactionID=rID;
    }
    public interface soslibbridge extends Library{
        soslibbridge bridge = (soslibbridge) Native.loadLibrary("/opt/lib/sosLibLinkv3.so", soslibbridge.class);
        
        ByValue runSolver3(String link,Pointer vals, int numVals, Pointer newparam, Pointer names,  int parametercount, Pointer rID);
        void outputCleanup3(Pointer pVals);
        
    }
    
    public double[] runsolver(double[] newparam) {
        soslibbridge runlib = soslibbridge.bridge;
        int num = varyingmetabolites.length+1;          //added extra ONE number to store the time steps needed to finish ode solution
        outputarray = new double[num];
        int vnamesize = 0;
        int pnamesize = 0;
        int rIDsize = 0;
        
        
        for(String name : varyingmetabolites){
            vnamesize+=name.length()+1;
        }
        
        for(int i=0; i<parametercount;i++){
            pnamesize+=parameterNames[i].length()+1;
            rIDsize+=reactionID[i].length()+1;
        }
        // Memory allocates a contiguous block of memory suitable for sharing with C directly - no copying is necessary
        Pointer parameter = new Memory(parametercount * Native.getNativeSize(Double.TYPE));
        Pointer variableNames = new Memory(vnamesize * Native.getNativeSize(Character.TYPE));
        Pointer parameterIDs = new Memory(pnamesize * Native.getNativeSize(Character.TYPE));
        Pointer reactionIDs = new Memory(rIDsize * Native.getNativeSize(Character.TYPE));
        // note: Memory instance (and its associated memory allocation) will be freed when ex8p goes out of scope
        for (int i=0; i<parametercount; i++) {
                parameter.setDouble(i * Native.getNativeSize(Double.TYPE), ((double)newparam[i]) );
        }
        long offset = 0;
        for (int i =0; i <varyingmetabolites.length; i++) {
            variableNames.setString(offset, varyingmetabolites[i]);                
            // null-terminator
            variableNames.setMemory(offset + varyingmetabolites[i].length(), 1, (byte)(0));
            offset += varyingmetabolites[i].length() + 1;
        }
       
        offset = 0;
        for (int i =0; i <parameterNames.length; i++) {
            parameterIDs.setString(offset, parameterNames[i]);
            // null-terminator
            parameterIDs.setMemory(offset + parameterNames[i].length(), 1, (byte)(0));
            offset += parameterNames[i].length() + 1;
        }
        offset = 0;
        for (int i =0; i <reactionID.length; i++) {
            reactionIDs.setString(offset, reactionID[i]);
            // null-terminator
            reactionIDs.setMemory(offset + reactionID[i].length(), 1, (byte)(0));
            offset += reactionID[i].length() + 1;
        }
        
        ByValue values = runlib.runSolver3(link, variableNames, num, parameter, parameterIDs, parametercount, reactionIDs);
        
        for (int i=0; i<num; i++) {
            // extract each double value from the buffer held by the struct
            outputarray[i] = values.vals.getDouble(i * Native.getNativeSize(Double.TYPE));
        }
        
        runlib.outputCleanup3(values.vals);
        values.clear();
        
        return outputarray;
    }
}
