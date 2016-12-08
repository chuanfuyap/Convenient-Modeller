/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Reaction;

import Species.*;
import java.util.ArrayList;
import org.sbml.jsbml.Reaction;

/**
 *
 * @author chuanfuyap
 */
public class MetabolicReaction {
    
    private ArrayList<Compound> substrates;
    private ArrayList<Compound> products;
    private ArrayList substratesStoichio;
    private ArrayList productsStoichio;
    private static int id = 0;
    private String ID;
    private int substrateStoichioIndex = 0;
    private int productStoichioIndex = 0;
    private String name;
    private Reaction reaction;
    private Enzyme enzyme = null;
    private ArrayList<Double> subconc;
    private ArrayList<Double> prodconc;
    private boolean proteinkinetics=false;
    private int regulation;
    private Compound modifier;
    
    public MetabolicReaction(String name){
        id = id + 1;
        this.name = name;
        this.ID = "R" + id;
        substrates = new ArrayList<Compound>();
        products = new ArrayList<Compound>();
        substratesStoichio = new ArrayList();
        productsStoichio = new ArrayList();
        subconc = new ArrayList<Double>();
        prodconc = new ArrayList<Double>();   
        regulation=1;
    }
    public MetabolicReaction (String name, int regulation, Compound modifier){
        this(name);
        this.regulation=regulation;
        this.modifier=modifier;
    }
    
    public MetabolicReaction(Reaction reaction){
        this.reaction = reaction;
        this.name = reaction.getName();
        ID = reaction.getId();
        id+=1;
        substrates = new ArrayList<Compound>();
        products = new ArrayList<Compound>();
        substratesStoichio = new ArrayList();
        productsStoichio = new ArrayList();
        subconc = new ArrayList<Double>();
        prodconc = new ArrayList<Double>();
        regulation=1;
        //1 is none, 2 is activator, 3 is inhibitor
    }
    
    public void setID(String ID){
        this.ID = ID;
    }
    
    public void addSubstrate(Compound compound) {
        substrates.add(compound);
        substratesStoichio.add(substrateStoichioIndex, 1);
        substrateStoichioIndex++;
        subconc.add(compound.getConcentration());
    }
    public void addSubstrates(ArrayList<Compound> substrates, ArrayList stoichio){
        for(int i =0; i <substrates.size();i++){
            addSubstrate(substrates.get(i),(int) stoichio.get(i));
        }
    }
    
    public void addSubstrate(Compound compound, int stoichio) {
        substrates.add(compound);
        substratesStoichio.add(substrateStoichioIndex, stoichio);
        substrateStoichioIndex++;
        subconc.add(compound.getConcentration());
    }
    
    public void addProduct(Compound compound) {
        products.add(compound);
        productsStoichio.add(productStoichioIndex, 1);
        productStoichioIndex++;
        prodconc.add(compound.getConcentration());
    }

    public void addProduct(Compound compound, int stoichio) {
        products.add(compound);
        productsStoichio.add(productStoichioIndex, stoichio);
        productStoichioIndex++;
        prodconc.add(compound.getConcentration());
    }
    public void addProducts(ArrayList<Compound> products, ArrayList stoichio){
        for(int i =0; i <products.size();i++){
            addProduct(products.get(i),(int) stoichio.get(i));
        }
    }
    public ArrayList getSubstrateStoichio() {
        return substratesStoichio;
    }
    
    public ArrayList<Compound> getSubstrates() {
        return substrates;
    }

    public ArrayList getProductStoichio() {
        return productsStoichio;
    }
        
    public ArrayList<Compound> getProducts() {
        return products;
    }
    
    public String getReactionID() {
        return ID;
    }
    
    public String getName() {
        return name;
    }
    
    public Reaction getReaction(){
        return reaction;
    }
    
    public Enzyme getReactionEnzyme(){
        return enzyme;
    }
    
    public void setReactionId(String id){
        ID = id;
    }
    
    public void setName(String n){
        name = n;
    }
    
    public void setEnzyme(Enzyme enz){
         enzyme = enz;
    }
    
    public void setPKinetics(boolean yesno){
        proteinkinetics=yesno;
    }
    
    public boolean getPKinetics(){
        return proteinkinetics;
    }
    
    public int getSubstrateStoicCo(int i){
        return Integer.parseInt(substratesStoichio.get(i).toString());
    }
    
    public int getProductStoicCo(int i){
       return Integer.parseInt(productsStoichio.get(i).toString()); 
    }
    public Compound getModifier(){
        return modifier;
    }
    
    public int getRegulation(){
        return regulation;
    }
    
    public void setRegulation(int num){
        regulation = num;
    }
    
    public void setModifier(Compound mod){
        modifier=mod;
    }
    
    public String toString() {
        String output = "";
        String subs = "";
        String prod = "";
        
        for (int i = 0; i < substrates.size(); i++) {
            subs += substratesStoichio.get(i) + substrates.get(i).getName();
            if (i!=substrates.size()-1){
                subs+=" + ";
            }
        }

        for (int i = 0; i < products.size(); i++) {
            prod += productsStoichio.get(i) + products.get(i).getName();
            if (i!=products.size()-1){
                prod += " + ";
            }
            
        }
        output += ("[" + subs + " <-----> " + prod + "]");
        
        return output;
    }
}
