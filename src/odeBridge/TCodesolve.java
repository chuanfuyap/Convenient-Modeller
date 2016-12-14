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
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

/**
 *
 * @author chuanfuyap
 */
public class TCodesolve {
    private String link;
    private int parametercount;
    private String[] varyingmetabolites;
    private String[] parameterNames;
    private String[] reactionID;
    private String[] reactionID2;
    private double[][] outputarray;
    
    
    public TCodesolve(String sbmllink, String[] variables, int count, String[] names, String[] rID, String[] rID2){
        link=sbmllink;
        parametercount=count;
        varyingmetabolites=variables;
        parameterNames=names;
        reactionID=rID;
        reactionID2=rID2;
    }
    public interface soslibbridge extends Library{
        soslibbridge bridge = (soslibbridge) Native.loadLibrary("/lib/SOSLib/dist/SOSlibJNIC.so", soslibbridge.class);
        
        void TCrunSolver(String link,Pointer vals, int numVals, Pointer newparam, Pointer names,  int parametercount, Pointer rID, int rIDsize, Pointer rID2, Pointer time, int last, PointerByReference RowArray, IntByReference RowCount);
        void TCarrayCleanup(Pointer RowArray, int count);
        
    }
    
    public double[][] runsolver(double[] newparam, double[] time) {
        soslibbridge runlib = soslibbridge.bridge;
        int num = varyingmetabolites.length+1;
        outputarray = new double[time.length][num+reactionID2.length];
        int vnamesize = 0;
        int pnamesize = 0;
        int rIDsize = 0;
        int rIDsize2 = 0;
        
        
        for(String name : varyingmetabolites){
            vnamesize+=name.length()+1;
        }
        for(int i=0; i<reactionID2.length;i++){
            rIDsize2+=reactionID2[i].length()+1;
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
        Pointer reactionIDs2 = new Memory(rIDsize2 * Native.getNativeSize(Character.TYPE));
        Pointer timepoints = new Memory(time.length * Native.getNativeSize(Double.TYPE));
        // note: Memory instance (and its associated memory allocation) will be freed when ex8p goes out of scope
        for (int i=0; i<parametercount; i++) {
            parameter.setDouble(i * Native.getNativeSize(Double.TYPE), ((double)newparam[i]) );
        }
        for(int i=0; i<time.length; i++){
            timepoints.setDouble(i * Native.getNativeSize(Double.TYPE), ((double)time[i]) );
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
        offset = 0;
        for (int i =0; i <reactionID2.length; i++) {
            reactionIDs2.setString(offset, reactionID2[i]);
            // null-terminator
            reactionIDs2.setMemory(offset + reactionID2[i].length(), 1, (byte)(0));
            offset += reactionID2[i].length() + 1;
        }
        
        PointerByReference rowPtrRef = new PointerByReference();
        IntByReference rowCountRef = new IntByReference();
        
        runlib.TCrunSolver(link, variableNames, num, parameter, parameterIDs, parametercount, reactionIDs, reactionID2.length, reactionIDs2, timepoints, time.length-1, rowPtrRef, rowCountRef);
        
        int rowcount = rowCountRef.getValue();
        Pointer rows = rowPtrRef.getValue();
        
        RowArray rowarray = new RowArray(rows);
        rowarray.read();
        
        RowArray[] rowarray2 = (RowArray[])rowarray.toArray(rowcount);
        int counter=num+reactionID2.length;
        for(int i=0;i<time.length;i++){
            for (int j=0; j<counter; j++) {
            //extract each double value from the buffer held by the struct
            outputarray[i][j] = rowarray2[i].vals.getDouble(j * Native.getNativeSize(Double.TYPE));
            }
        }
        
        runlib.TCarrayCleanup(rows, time.length);
                
        return outputarray;
    }
        
}
