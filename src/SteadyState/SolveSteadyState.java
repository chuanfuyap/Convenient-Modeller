package SteadyState;

import Main.ModelReaction;
import Main.SBMLFile;
import Species.Compound;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import odeBridge.odesolveSS;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author chuanfuyap
 */
public class SolveSteadyState {
    
    private double Tolerance = 5E-7;
    private double NRDTolerance1 = 5E-3;
    private double NRDTolerance2 = 5E-5;
    private int maxiteration = 100;
    private int dampingiteration = 200;
    private ArrayList<Compound> Compounds;
    private ArrayList<ModelReaction> Reactions;
    private double[] solvedvector;
    private double[] Parameters;
    private double[] dampedvector;
    private double[] dampedvector2;
    private double[] odesolvedvalue;
    private HashMap dampedMap = new HashMap<>();
    private HashMap solvedMap = new HashMap<>();
    private boolean secondrun = false;
    private String link;
    private boolean solution=false;
    
    //the SSvalue input here is for when doing parameter estimation it would be close to the ideal solution
    public SolveSteadyState(ArrayList<Compound> Compounds, ArrayList<ModelReaction> Reactions, double[] Parameters, String link) {
        this.Compounds=Compounds;
        this.Reactions=Reactions;
        this.Parameters=Parameters;                             //parameters being solved
        this.link=link;
    }
    
    public boolean solveSS() throws InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, ModelOverdeterminedException, XMLStreamException, IOException {
        File SBMLFileLink = new File(link);                     //needed to construct odesolver object
        SBMLFile sbmlfile = new SBMLFile(SBMLFileLink);         //needed to construct odesolver object
        String[] compID = ActiveCompoundIDs();
        String[] paranames = sbmlfile.getParameterNames();
        String[] reactionID = sbmlfile.getReactionID();
        
        odesolveSS ode = new odesolveSS(link, compID, Parameters.length, paranames, reactionID);        //things needed to solve steady state
        
        odesolvedvalue = ode.runsolver(Parameters);
                
        VariableVector VV = new VariableVector(Compounds, odesolvedvalue);                              //things needed to solve steady state
        
        //RUN NEWTONRAPHSON WITH DAMPING, it will call odesolver up to 10000s if it fails to bring it towards a solution
        boolean dampsolution=false;
        boolean NRDSolve=false;
        //runs NR with damping
        try{
            NRDSolve=NewtonRaphsonWithDamping(VV, NRDTolerance1);
        }catch(Exception e){
            if(dampedvector==null){
                dampedvector = new double[compID.length];
                for (int k = 0 ; k<compID.length; k++){
                    dampedvector[k]=0;
                }
            }
        }
//        System.out.println(NRDSolve);
        if (NRDSolve==true){
            dampsolution=true;
        }
        
        //AFTER DAMPING is done, it continues with normal newton raphson to refine answer
        if (dampsolution==false){
            //in the event that damping failed, it would be null, thus nothing happens, the hashmaps would be empty.
            solution=false;
        }else{
            solution=true;
            //this happens when damping is SUCCESSFUL
            VariableVector newvector = new VariableVector(Compounds, dampedvector);            
            //trying normal newton raphson
            boolean NRsolve = false;
            try{
                NRsolve=NewtonRaphson(newvector);
            }catch(Exception e){
                VariableVector newervector = new VariableVector(Compounds, dampedvector);
                try{
                    secondrun=true;
                    NRsolve=NewtonRaphsonWithDamping(newervector, NRDTolerance2);
                }catch(Exception e2){
                }
                
            }
//            System.out.println(NRsolve);
            if(NRsolve==false){
                solvedvector=dampedvector;
                solvedMap=dampedMap;
            }else{                
                for (int j =0; j<solvedvector.length; j++){
                    solvedMap.put(compID[j], solvedvector[j]);
                }
            }
        }
        return solution;
    }
    
    //NEWTONRAPHSON WITH INCREASING DAMPING METHOD
    private boolean NewtonRaphsonWithDamping(VariableVector Variables, double DampingTolerance){
        boolean canitbesolved=false;
        boolean solved=false;
        boolean NaN=false;
        
        FunctionVector Functions = new FunctionVector(Variables, Compounds, Reactions, Parameters);
        JacobianMatrix Jacobian = new JacobianMatrix(Variables, Compounds, Reactions, Parameters);
        
        RealMatrix vvector = Variables.getVector();
        RealMatrix fvector = Functions.getVector();        
        RealMatrix inverse = MatrixUtils.inverse(Jacobian.getMatrix());
        
        double dampvalue=0;
        for (double i = 0; i<dampingiteration && solved==false && NaN==false; i++){
            if(dampvalue<32.0){
                dampvalue = i/4;
            }
            double damping = Math.pow(2,-dampvalue);
            RealMatrix damp = inverse.multiply(fvector).scalarMultiply(damping);
            RealMatrix x1 = vvector.subtract(damp);
            double[] newvalue= x1.getColumn(0);
            
            for (double num : newvalue){                //in case theres error which results in NULL values, which would return the method has failed
                NaN=Double.isNaN(num);
            }
            
            VariableVector VV2 = new VariableVector(Compounds, newvalue);
            FunctionVector FV2 = new FunctionVector(VV2, Compounds, Reactions, Parameters);
            JacobianMatrix JM2 = new JacobianMatrix(VV2, Compounds, Reactions, Parameters);

            vvector = VV2.getVector();
            fvector = FV2.getVector();
            inverse = MatrixUtils.inverse(JM2.getMatrix());

            if (secondrun==false) {
                if (IsNRDSolved(FV2, VV2, DampingTolerance) == true) {
                    solved = true;
                    canitbesolved = true;
                    dampedvector = vvector.getColumn(0);
                    int count = 0;
                    for (Compound comp : Compounds) {
                        if (comp.getBoundaryCondition() == false) {
                            dampedMap.put(comp.getID(), dampedvector[count]);
                            count++;
                        }
                    }
//                    System.out.println("damping done");
                }
            }else{
                String[] compID = ActiveCompoundIDs();
                if (IsNRDSolved(FV2, VV2, DampingTolerance) == true) {
                    solved=true;
                    canitbesolved = true;
                    dampedvector2=vvector.getColumn(0);
                    solvedvector=dampedvector2;
                    for (int j =0; j<dampedvector2.length; j++){
                        solvedMap.put(compID[j], dampedvector2[j]);
                    }
                }
            }
        }
        return canitbesolved;
    }
    
    //NORMAL NEWTONRAPHSON METHOD
    private boolean NewtonRaphson(VariableVector Variables){
        boolean canitbesolved=false;
        boolean solved=false;
        boolean NaN=false;
        FunctionVector Functions = new FunctionVector(Variables, Compounds, Reactions, Parameters);
        JacobianMatrix Jacobian = new JacobianMatrix(Variables, Compounds, Reactions, Parameters);
        
        RealMatrix vvector = Variables.getVector();
        RealMatrix fvector = Functions.getVector();        
        RealMatrix inverse = MatrixUtils.inverse(Jacobian.getMatrix());
        
        for (int i = 0; i<maxiteration && solved==false && NaN==false; i++){
            RealMatrix x1 = vvector.subtract(inverse.multiply(fvector));
            double[] newvalue= x1.getColumn(0);
            
            for (double num : newvalue){                            //in case theres error which results in NULL values, which would return the method has failed
                NaN=Double.isNaN(num);
            }
            
            VariableVector testvvector = new VariableVector (Compounds, newvalue);
            FunctionVector testfvector = new FunctionVector(testvvector, Compounds, Reactions, Parameters);
            solved = NRSolve(testfvector, newvalue);
            
            VariableVector VV2 = new VariableVector (Compounds, newvalue);
            FunctionVector FV2 = new FunctionVector(VV2, Compounds, Reactions, Parameters);
            JacobianMatrix JM2 = new JacobianMatrix(VV2, Compounds, Reactions, Parameters);
            vvector = VV2.getVector();
            fvector = FV2.getVector();
            inverse = MatrixUtils.inverse(JM2.getMatrix());
            
            if(solved==true){
                canitbesolved=true;
            }
            
        }
        return canitbesolved;
    }
    
    //METHOD to store all the IDs of compounds that are active to be used when extracting variable conc from ODEsolver
    private String[] ActiveCompoundIDs(){
        
        int numComp = 0;
        for (Compound name : Compounds){
            if (name.getBoundaryCondition()==false){
                numComp++;
            }
        }
        
        String[] activecompID = new String[numComp]; //species ID that are varying
        int tick = 0;
        for (Compound name : Compounds){
            if (name.getBoundaryCondition()==false){
                activecompID[tick]=name.getID();
                tick++;
            }
        }
        return activecompID;
    }
    
    //check if the NewtonRaphson method solution has worked
    private boolean IsNRDSolved(FunctionVector FunctionVector, VariableVector VVector, double Tolerance){
        
        boolean solved=false;
        
        RealMatrix tempmatrix = FunctionVector.getVector();
        
        double[] values = tempmatrix.getColumn(0);
//        System.out.println("damping");
        double total=0;
        for (double num : values){
            total+=Math.abs(num*num);
        }
        
        double avg = total/values.length;
//        System.out.println(avg+"\n");
        if (Math.abs(avg) <= Tolerance){
            solved=true;
        }
        
        return solved;
    }
    
    private boolean NRSolve(FunctionVector FunctionVector, double[] variablevector){
        RealMatrix tempmatrix = FunctionVector.getVector();
        double[] functionvalue = tempmatrix.getColumn(0);
        boolean solved=false;
        
        double total=0;
        for(double num : functionvalue){
            total+=Math.abs(num*num);
        }
        double avg = total/functionvalue.length;
//        System.out.println(avg+"\n");
        if(avg<=Tolerance){
            solved=true;
            solvedvector=variablevector;
        }
        
        return solved;
    }
    
    public double[] getSolvedVector(){
        return solvedvector;
    }
    
    public HashMap getSolvedMMap(){
        return solvedMap;
    }
    
    public HashMap getSolvedFMap(){
        
        VariableVector vv = new VariableVector(Compounds, solvedvector);
        FunctionBuilder fb = new FunctionBuilder(vv, Reactions, Parameters);
        
        return fb.getSolvedFunctions();
        
    }
}
