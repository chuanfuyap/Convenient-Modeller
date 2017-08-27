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

/**
 *
 * @author chuanfuyap
 */
public class FunctionBuilder {
    
    private HashMap AllFunctions = new HashMap<>();
    private HashMap solvedFunctions = new HashMap<>();
    
    //builds function for all the reactions
    public FunctionBuilder (VariableVector VVector, ArrayList<ModelReaction> Reactions, double[] Parameters){
        
        HashMap metmap = VVector.getMetMap();
        int Tracker = 0;
        
        for (ModelReaction Reaction : Reactions){
            int numSub = Reaction.getSubstrates().size();
            int numProd = Reaction.getProducts().size();
            int numPara=2;
            
            int regulation = Reaction.getRegulation();
            
            if(regulation==2||regulation==3){
                numPara+=1;
            }else if (regulation == 4){
                numPara+=2;
            }
            numPara+=numSub;
            numPara+=numProd;
            
            CKFunction Function;
            ArrayList<Compound> subs = Reaction.getSubstrates();
            ArrayList<Compound> prod = Reaction.getProducts();
            ArrayList<Compound> mods = Reaction.getModifier();
            
            String ReactionName = Reaction.getName();
            String ReactionID = Reaction.getReactionID();
            
            ArrayList<Double> subconc = new ArrayList<>();
            ArrayList<Double> prodconc = new ArrayList<>();
            ArrayList<Double> modconc = new ArrayList<>();
            ArrayList<Double> kpvalue = new ArrayList<>();
            
            double enzconc = Reaction.getEnzyme().getConcentration();
            
            for (Compound substrate : subs){
                subconc.add((Double) metmap.get(substrate.getName()));
            }
            
            for (Compound product : prod){
                prodconc.add((Double) metmap.get(product.getName()));
            }
            if(regulation>1){
            for (Compound modifier : mods){
                modconc.add((Double) metmap.get(modifier.getName()));
            }
            }
            
            for (int i=0; i < numPara; i++){
                kpvalue.add(Parameters[Tracker]);
                Tracker++;
            }
            
            ArrayList substoichio = Reaction.getSubStoichio();
            ArrayList prodstoichio = Reaction.getProdStoichio();
            
//            System.out.println(Reaction.getName());
            Function = new CKFunction(subconc, prodconc, kpvalue, enzconc, regulation, modconc, substoichio, prodstoichio);
            
            AllFunctions.put(ReactionName, Function);
            solvedFunctions.put(ReactionID, Function.SolveFunction());
//            System.out.println(ReactionID+"\t"+Function.SolveFunction());
        }
        
    }
    
    public HashMap returnFunctions(){
        return AllFunctions;
    }
    
    public HashMap getSolvedFunctions(){
        return solvedFunctions;
    }

}
