/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SteadyState;

import Main.ModelReaction;
import Species.Compound;
import java.util.ArrayList;

/**
 *
 * @author chuanfuyap
 */
public class MetaboliteRateFunction {
    
    private ArrayList stoichio1 = new ArrayList<>();
    private ArrayList stoichio2 = new ArrayList<>();
    private ArrayList<String> reactionname1 = new ArrayList<>(); 
    private ArrayList<String> reactionname2 = new ArrayList<>(); 
    private ArrayList<ArrayList> reaction1sub= new ArrayList<>(); ;
    private ArrayList<ArrayList> reaction1prod= new ArrayList<>(); ;
    private ArrayList<ArrayList> reaction2sub= new ArrayList<>(); ;
    private ArrayList<ArrayList> reaction2prod= new ArrayList<>(); ;
    private String name;
    private ArrayList<String> allcompounds = new ArrayList<>();
    private ArrayList<String> reaction1compounds = new ArrayList<>();
    private ArrayList<String> reaction2compounds = new ArrayList<>();
    private ArrayList<String> reactionname1id = new ArrayList<>(); 
    private ArrayList<String> reactionname2id = new ArrayList<>(); 
    private ArrayList<String> modifiers = new ArrayList<>(); 
    
    //this class is to build each function for a given compound. 
    public MetaboliteRateFunction(Compound Compound, ArrayList<ModelReaction> Reactions){
        
        //function to be built for this compound, eg: d[species]/dt=production-consumption
        name = Compound.getName();
                
        //reads into every reaction to look for which reaction the compound is participating in. 
        for (ModelReaction Reaction : Reactions){
            
            ArrayList<Compound> subs = Reaction.getSubstrates();
            ArrayList<Compound> prod = Reaction.getProducts();
            ArrayList<Integer> substoichio = Reaction.getSubstratesStoichio();
            ArrayList<Integer> prodstoichio = Reaction.getProductsStoichio();
            ArrayList<String> subnames =new ArrayList<>();
            ArrayList<String> prodnames =new ArrayList<>();           
            if(Reaction.getRegulation()>1){
                if(Reaction.getRegulation()<4){
                    modifiers.add(Reaction.getModifier().get(0).getName());
                    allcompounds.add(Reaction.getModifier().get(0).getName());
                }else if (Reaction.getRegulation()==4){
                    modifiers.add(Reaction.getModifier().get(0).getName());
                    allcompounds.add(Reaction.getModifier().get(0).getName());
                    modifiers.add(Reaction.getModifier().get(1).getName());
                    allcompounds.add(Reaction.getModifier().get(1).getName());
                }
            }
            
            
            for (Compound sub : subs){
                subnames.add(sub.getName());
            }
            
            for (Compound product : prod){
                prodnames.add(product.getName());
            }
            
            //if it is a product in this reaction, it means production, thus POSITIVE/reaction1
            for (int i = 0; i < prod.size(); i++){
                String tempname = prod.get(i).getName();
                if (name.equals(tempname)){
                    
                    stoichio1.add(prodstoichio.get(i));
                    reactionname1.add(Reaction.getName());
                    reaction1sub.add(subnames);
                    reaction1prod.add(prodnames);
                    allcompounds.addAll(subnames);
                    allcompounds.addAll(prodnames);
                    reaction1compounds.addAll(subnames);
                    reaction1compounds.addAll(prodnames);
                    reactionname1id.add(Reaction.getReactionID());
                }
            }
            
            //if it is a substrate in this reaction, it means consumption, thus NEGATIVE/reaction2
            for (int j = 0; j < subs.size(); j++){
                String tempname = subs.get(j).getName();
                if (name.equals(tempname)){
                    stoichio2.add(substoichio.get(j));
                    reactionname2.add(Reaction.getName());
                    reaction2sub.add(subnames);
                    reaction2prod.add(prodnames);
                    allcompounds.addAll(subnames);
                    allcompounds.addAll(prodnames);
                    reaction2compounds.addAll(subnames);
                    reaction2compounds.addAll(prodnames);
                    reactionname2id.add(Reaction.getReactionID());
                } 

            }
            
        }
    }
    
    public void returnFunction(){
        System.out.println(stoichio1+"*"+reactionname1+"\t-\t"+stoichio2+"*"+reactionname2);
    }
    
    public ArrayList<String> getProduction(){
        return reactionname1;
    }
    
    public ArrayList<String> getConsumption(){
        return reactionname2;
    }
    
    public ArrayList prodstoichio(){
        return stoichio1;
    }
    
    public ArrayList consumpstoichio(){
        return stoichio2;
    }
    
    public ArrayList<ArrayList> getreaction1sub(){
        return reaction1sub;
    }
    
    public ArrayList<ArrayList> getreaction1prod(){
        return reaction1prod;
    }
    
    public ArrayList<ArrayList> getreaction2sub(){
        return reaction2sub;
    }
    
    public ArrayList<ArrayList> getreaction2prod(){
        return reaction2prod;
    }
    
    public String getName(){
        return name;
    }
    
    public ArrayList<String> getAllCompounds(){
        return allcompounds;
    }
    
    public ArrayList<String> getProdCompounds(){
        return reaction1compounds;
    }
    
    public ArrayList<String> getConsCompounds(){
        return reaction2compounds;
    }
    
    public ArrayList<String> getProductionID(){
        return reactionname1id;
    }
    
    public ArrayList<String> getConsumptionID(){
        return reactionname2id;
    }
    
    public ArrayList<String> getmodifiers(){
        return modifiers;
    }
}
