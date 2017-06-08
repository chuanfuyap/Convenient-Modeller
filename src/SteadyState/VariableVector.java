/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SteadyState;

import Species.Compound;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 *
 * @author chuanfuyap
 */
public class VariableVector {
    
    int numMet=0, numInOut=0; //comp = active metabolite, InOut=input/output where boundary condition is true.
    HashMap metMap = new HashMap<String, Double>();
    String[] metName, inoutName;
    double[] metConc, inoutConc;
    RealMatrix VVector; //Variable Vector
    
    public VariableVector(ArrayList<Compound> Compounds, double[] activeConc){
        
        for (Compound Compound : Compounds) {
            if (Compound.getBoundaryCondition()==false){
                numMet++;
                
            }
        }
        
        for (Compound Compound : Compounds) {
            if (Compound.getBoundaryCondition()==true){
                numInOut++;
            }
        }
        
        metName = new String[numMet];
        inoutName = new String[numInOut];
        metConc = new double[numMet];
        inoutConc = new double[numInOut];
        int counter =0;
        for (int i=0; i<Compounds.size(); i++) {
            
            if (Compounds.get(i).getBoundaryCondition()==false ){
                metName[counter]=Compounds.get(i).getName();
                metConc[counter]=activeConc[counter];
//                System.out.println(metConc[counter]);
                metMap.put(Compounds.get(i).getName(), activeConc[counter]);
                counter++;
                
            }
        } 
        
        counter =0;
        for (int i=0; i<Compounds.size(); i++) {
            if (Compounds.get(i).getBoundaryCondition()==true ){
                inoutName[counter]=Compounds.get(i).getName();
                if(Compounds.get(i).getBoundaryCondition()==true){
                    inoutConc[counter]=Compounds.get(i).getConcentration();
                    metMap.put(Compounds.get(i).getName(), Compounds.get(i).getConcentration());
                }
                counter++;
            }
        }
        try{
            VVector = MatrixUtils.createColumnRealMatrix(metConc);
        }catch(Exception e){
            System.out.println(e);
            for(double eh : metConc){
                System.out.println(eh);
            }
            
        }
        
    }
    
        
    public String[] getMetName(){
        return metName;
    }
    
    public String[] getInOutName(){
        return inoutName;
    }
    
    public double[] getMetConc(){
        return metConc;
    }
    
    public double[] getInOutConc(){
        return inoutConc;
    }
    
    public HashMap getMetMap(){
        return metMap;
    }
    
    public RealMatrix getVector(){
        return VVector;
    }
}
