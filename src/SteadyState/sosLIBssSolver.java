package SteadyState;

import Main.ModelReaction;
import Main.SBMLFile;
import Species.Compound;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import odeBridge.SSodesolve;
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
public class sosLIBssSolver {
    
    private ArrayList<Compound> Compounds;
    private ArrayList<ModelReaction> Reactions;
    private double[] solvedvector;
    private double[] Parameters;
    private HashMap solvedMap = new HashMap<>();
    private String link;
    private boolean solution=false;
    private double[] initialValues;
    
    //the SSvalue input here is for when doing parameter estimation it would be close to the ideal solution
    public sosLIBssSolver(ArrayList<Compound> Compounds, ArrayList<ModelReaction> Reactions, double[] Parameters, String link) {
        this.Compounds=Compounds;
        this.Reactions=Reactions;
        this.Parameters=Parameters;                             //parameters being solved
        this.link=link;
        this.initialValues=getInitialConc();
    }
    
    public boolean solveSS() throws InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, ModelOverdeterminedException, XMLStreamException, IOException {
        File SBMLFileLink = new File(link);                     //needed to construct odesolver object
        SBMLFile sbmlfile = new SBMLFile(SBMLFileLink);         //needed to construct odesolver object
        String[] compID = ActiveCompoundIDs();
        String[] paranames = sbmlfile.getParameterNames();
        String[] reactionID = sbmlfile.getReactionID();
        
        SSodesolve ode = new SSodesolve(link, compID, Parameters.length, paranames, reactionID);        //things needed to solve steady state
        
        solvedvector = ode.runsolver(Parameters);
        
        double total=0;
        for(int i=0;i<solvedvector.length-1;i++){
            total+=solvedvector[i];
        }
        if(total>0){//if soslib has error total would be 0 as theres no output, resulting in false solution, therefore not needing to move forward
        double timepointNeeded = solvedvector[solvedvector.length-1]; //solvedvector includes an extra index which is the number of timepoints taken to solve the model. 
        
        if(timepointNeeded<1001.0){
            double[]tempvector = new double[solvedvector.length-1];
            for (int j =0; j<solvedvector.length-1; j++){
                solvedMap.put(compID[j], solvedvector[j]);
                tempvector[j]=solvedvector[j];
            }
            solvedvector=tempvector;
            solution=wasitGood();       //added check for when there is CVODE error which gives me back initial values and a timepointNeed of 1
        }
        }
        return solution;
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
    
    private double[] getInitialConc(){
        int numComp = 0;
        for (Compound name : Compounds){
            if (name.getBoundaryCondition()==false){
                numComp++;
            }
        }
        initialValues = new double[numComp];
        int tick = 0;
        for (Compound comp : Compounds){
            if (comp.getBoundaryCondition()==false){
                initialValues[tick]=comp.getConcentration();
                tick++;
            }
        }
        return initialValues;
    }
    
    private boolean wasitGood(){
        boolean good = true;
        double total=0;
        for(int i =0; i<initialValues.length;i++){
            double difference=initialValues[i]-solvedvector[i];
            total+=difference;
        }
        
        if(total==0){
            good=false;
        }
        
        return good;
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
