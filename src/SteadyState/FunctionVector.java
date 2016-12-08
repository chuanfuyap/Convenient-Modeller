/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SteadyState;

import Main.ModelReaction;
import Species.Compound;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 *
 * @author chuanfuyap
 */
public class FunctionVector {
    
    private RealMatrix FVector;
    private double[] Fvalues;    
    int numComp;
    HashMap FuncMap = new HashMap<>();
    
    //this class is to build functions for all the boundarycondition=false compound and to put them into a vector.
    public FunctionVector (VariableVector VVector, ArrayList<Compound> Compounds, ArrayList<ModelReaction> Reactions, double[] Parameters){
        
        numComp = Compounds.size();
        
        FunctionBuilder allfunctions = new FunctionBuilder(VVector, Reactions, Parameters);
        double[] apple = VVector.getMetConc();
                
        HashMap functionmap = allfunctions.returnFunctions();
        
        for (Compound Compound : Compounds) {
            
            if (Compound.getBoundaryCondition()==true ){
                numComp--;
            }
        
        }
        
        Fvalues = new double[numComp];
        
        int counter=0;
        for (int i = 0; i < Compounds.size(); i++) {
            
            if (Compounds.get(i).getBoundaryCondition()==false ){
                
                MetaboliteRateFunction MRF = new MetaboliteRateFunction(Compounds.get(i), Reactions);
                
                ArrayList<CKFunction> tempfunc1 = new ArrayList<>();
                ArrayList<CKFunction> tempfunc2 = new ArrayList<>();
                ArrayList prodstoichio = MRF.prodstoichio();
                ArrayList consumpstoichio = MRF.consumpstoichio();
                for (String production : MRF.getProduction()) {
                    tempfunc1.add((CKFunction) functionmap.get(production));
                }
                for (String consumption : MRF.getConsumption()) {
                    tempfunc2.add((CKFunction) functionmap.get(consumption));                    
                }
                double production=0;
                for(int j = 0 ; j < tempfunc1.size();j++){
                    production+=((int) prodstoichio.get(j))*tempfunc1.get(j).SolveFunction();         
                }
                double consumption=0;
                for(int j = 0 ; j < tempfunc2.size();j++){
                    consumption+=((int) consumpstoichio.get(j))*tempfunc2.get(j).SolveFunction();
                }
                
                double value = production-consumption;
                Fvalues[counter] = value;
//                System.out.println(production +"\t"+consumption+"=\t"+value);
                FuncMap.put(Compounds.get(i).getName(),value);
                
                counter++;
                
            }
            
        }
        
        FVector = MatrixUtils.createColumnRealMatrix(Fvalues);
        
    }
    
    public RealMatrix getVector(){
        return FVector;
    }
    
    public double[] getVectorArray(){
        return Fvalues;
    }
    
    public HashMap getFuncMap(){
        return FuncMap;
    }
}
