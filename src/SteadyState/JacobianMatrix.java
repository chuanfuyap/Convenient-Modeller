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
public class JacobianMatrix {
    
    private int numComp;
    private double[][] matrix;  //matrix
    private String[][] stringmatrix;
    private ArrayList<MetaboliteRateFunction> MRFlist = new ArrayList<>();  //stores all metabolite function (production-consumption)
    private ArrayList<String> compoundname = new ArrayList<>();  //stores compoundname
    private HashMap MRFmap = new HashMap<>();     //maps compound to their respective function
    private RealMatrix Jmatrix;
    private ArrayList<ModelReaction> Reactions;
    
    public JacobianMatrix (VariableVector Vvector, ArrayList<Compound> Compounds, ArrayList<ModelReaction> Reactions, double[] Parameters){
        
        this.Reactions = Reactions;
        FunctionBuilder allfunctions = new FunctionBuilder(Vvector, Reactions, Parameters);
        HashMap functionmap = allfunctions.returnFunctions();   //build the Convenience Kinetics Function for all the enzymatic reactions
        HashMap mod2reactionmap = modifier2reaction(Reactions);
        
        numComp = Compounds.size();
        
        //determine number compoounds that are a variable boundarycondition=falses
        for (Compound Compound : Compounds) {
            if (Compound.getBoundaryCondition()==true){
                numComp--;
            }
        }
        
        //instantiate the matrix
        matrix = new double[numComp][numComp];
        stringmatrix = new String[numComp][numComp];
        
        //generates all the metabolite rate function and maps them to their name
        for (int i = 0; i < Compounds.size(); i++) {
            if (Compounds.get(i).getBoundaryCondition()==false){
                MetaboliteRateFunction MRF = new MetaboliteRateFunction(Compounds.get(i), Reactions);
                compoundname.add(Compounds.get(i).getName());
                MRFmap.put(Compounds.get(i).getName(), MRF);
            }
        }
        
        for (int row =0; row < compoundname.size(); row++){                     //row
            String rowname=compoundname.get(row);   //rowname here would be the metabolites rate of change function name
            for (int col =0; col <compoundname.size(); col++ ){                 //column                
                String colname=compoundname.get(col);
                MetaboliteRateFunction rowMRF = (MetaboliteRateFunction) MRFmap.get(rowname);
//                ArrayList<CKFunction> productionFunction = new ArrayList<>();
//                ArrayList<CKFunction> consumptionFunction = new ArrayList<>();
                ArrayList stoichio1 = rowMRF.prodstoichio();
                ArrayList stoichio2 = rowMRF.consumpstoichio();
                
                ArrayList<ArrayList> r1sub,r1prod,r2sub,r2prod;
                r1sub = rowMRF.getreaction1sub();
                r1prod = rowMRF.getreaction1prod();
                r2sub = rowMRF.getreaction2sub();
                r2prod = rowMRF.getreaction2prod();
                ArrayList<String> mods = rowMRF.getmodifiers();
                
                ArrayList<String> reaction1comp = rowMRF.getProdCompounds();
                ArrayList<String> reaction2comp = rowMRF.getConsCompounds();
                ArrayList<String> allcompound = rowMRF.getAllCompounds();
                
                double overall=0; //the overall value of this element in matrix, 0 if variable not involved in row's function.
                double prodreaction=0;
                double consumereaction=0;
                String overalls="";
                String prodreactions="";
                String consumereactions="";
                
                //loop to determine if this column's variable is in this row's function
                boolean isinRowFunction=false;
                boolean mod = false;
                for (String name : allcompound){
                    if (name.equals(colname)){
                        isinRowFunction=true;
                    }                    
                }
                
                for (String name : mods){
                    if (name.equals(colname)){
                        mod=true;
                    }                    
                }
                
                if (isinRowFunction==false){
                    overall+=0;
                }else{  //if it is in the Row, do, search to see if the production reaction(R1) or consumption reaction(R2)
                         
                    boolean inR1 = false;
                    boolean inR2 = false;

                    for (String r1name : reaction1comp){                        
                        if (r1name.equals(colname)){
                           inR1=true;
                        }
                    }
                    for (String r2name : reaction2comp){
                        if (r2name.equals(colname)){
                            inR2=true;
                        }
                    }
                    
                    if(mod == false ){     //if column name is not a modifier
                        
                        if(inR1==false){    //if it is not involved in production it becomes 0-consumption
                            prodreaction+=0;
                            prodreactions+=0;
                        }else{              //now search whether its the substrate or product in the production reaction 

                            for(int j=0; j < r1sub.size(); j++ ){
                                for(int i =0; i < r1sub.get(j).size(); i ++){       //loop to look through all the reactions involved production side of the rate function
                                    if(r1sub.get(j).get(i).equals(colname)){        //loop through all the SUBSTRATES involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(j));                                    
                                        prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivate(true, i);
                                        prodreactions+=(int) stoichio1.get(j)+"Ps"+(i+1);
                                    }
                                }
                            }
                            for(int j=0; j < r1prod.size(); j++ ){
                                for(int i =0; i < r1prod.get(j).size(); i ++){      //loop to look through all the reactions involved production side of the rate function
                                    if(r1prod.get(j).get(i).equals(colname)){       //loop through all the PRODUCTS involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(j));
                                        prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivate(false, i);
                                        prodreactions+=(int) stoichio1.get(j)+"Pp"+(i+1);
                                    }
                                }
                            }
                        }
                        
                        if(inR2==false){    //if it is not involved in consumption reaction (production - 0)
                            consumereaction+=0;
                            consumereactions+=0;
                        }else{
                            for(int j=0; j < r2sub.size(); j++ ){
                                for(int i =0; i < r2sub.get(j).size(); i ++){       //loop to look through all the reactions involved consumption side of the rate function
                                    if(r2sub.get(j).get(i).equals(colname)){        //loop through all the SUBSTRATES involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(j));
                                        consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivate(true, i);
                                        consumereactions+=(int) stoichio2.get(j) +"Cs"+(i+1);
                                    }
                                }
                            }
                            for(int j=0; j < r2prod.size(); j++ ){
                                for(int i =0; i < r2prod.get(j).size(); i ++){      //loop to look through all the reactions involved consumption side of the rate function
                                    if(r2prod.get(j).get(i).equals(colname)){       //loop through all the PRODUCTS involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(j));
                                        consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivate(false, i);
                                        consumereactions+=(int) stoichio2.get(j) +"Cp"+(i+1);
                                    }
                                }
                            }
                        }
                        
                    }else{      //if column name is a modifier
                        
                        boolean inR1m = false;
                        boolean inR2m = false;
                        
                        for(int i =0;i<rowMRF.getProduction().size();i++){
                            String reactionName = rowMRF.getProduction().get(i);
                            ArrayList<Compound> rowModifier = (ArrayList) mod2reactionmap.get(reactionName);
                            for(Compound rowMod : rowModifier){
                                if(rowMod.getName().equals(colname)){
                                    inR1m=true;
                                }
                            }
                        }
                        
                        for(int i =0;i<rowMRF.getConsumption().size();i++){
                            String reactionName = rowMRF.getConsumption().get(i);
                            ArrayList<Compound> rowModifier = (ArrayList) mod2reactionmap.get(reactionName);
                            for(Compound rowMod : rowModifier){
                                if(rowMod.getName().equals(colname)){
                                    inR2m=true;
                                }
                            }
                        }
                        
                        if(inR1m==false&&inR1==false){    //if it is not involved in production it becomes 0-consumption
                            prodreaction+=0;
                            prodreactions+=0;
                        }  //now search whether its a sub/prod/mod/mod+sub/mod+prod in the production reaction 
                        
                        else if(inR1m==true&&inR1==false){      // just mod
                            
                            for(int i =0;i<rowMRF.getProduction().size();i++){
                                if(isModInCol(mod2reactionmap, rowMRF.getProduction().get(i), colname)){
                                    CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(i));
                                    if(isitDoubleModifier(rowMRF.getProduction().get(i))==false){
                                        prodreaction+= (int) stoichio1.get(i) * function.getPartialDerivateWithMod();
                                        prodreactions+="justmod";
                                    }else{
                                        prodreaction+= (int) stoichio1.get(i)* function.getPartialDerivateWithDoubleMod(isitActivatorOrInhibitor(rowMRF.getProduction().get(i), colname));
                                        prodreactions+="justDmod";
                                    }
                                }
                            }
                        }
                        
                        else if(inR1m==false&&inR1==true){     // just sub/prod
                            for(int j=0; j < r1sub.size(); j++ ){
                                for(int i =0; i < r1sub.get(j).size(); i ++){       //loop to look through all the reactions involved production side of the rate function
                                    if(r1sub.get(j).get(i).equals(colname)){        //loop through all the SUBSTRATES involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(j));
                                        prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivate(true, i);
                                        prodreactions+=(int) stoichio1.get(j)+"Ps"+(i+1);
                                    }
                                }
                            }
                            for(int j=0; j < r1prod.size(); j++ ){
                                for(int i =0; i < r1prod.get(j).size(); i ++){      //loop to look through all the reactions involved production side of the rate function
                                    if(r1prod.get(j).get(i).equals(colname)){       //loop through all the PRODUCTS involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(j));
                                        prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivate(false, i);
                                        prodreactions+=(int) stoichio1.get(j)+"Pp"+(i+1);
                                    }
                                }
                            }
                        }
                        
                        else if(inR1m==true&&inR1==true){      // mod + sub/prod
                            for(int j=0; j < r1sub.size(); j++ ){
                                for(int i =0; i < r1sub.get(j).size(); i ++){       //loop to look through all the reactions involved production side of the rate function
                                    if(r1sub.get(j).get(i).equals(colname)){        //loop through all the SUBSTRATES involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(j));
                                        
                                        if(isModInCol(mod2reactionmap, rowMRF.getProduction().get(j), colname)){
                                            if(isitDoubleModifier(rowMRF.getProduction().get(j))==false){
                                            prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivateWithModSubProd(true, i);
                                            prodreactions+=(int) stoichio1.get(j)+"Pms"+(i+1);
                                            }else{
                                                prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivateWithDoubleModSubProd(isitActivatorOrInhibitor(rowMRF.getProduction().get(j), colname), true, i);
                                                prodreactions+=(int) stoichio1.get(j)+"PDms"+(i+1);
                                            }
                                        }else{
                                            prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivate(true, i);
                                            prodreactions+=(int) stoichio1.get(j)+"Ps"+(i+1);
                                        }                                        
                                    }
                                }
                            }
                            for(int j=0; j < r1prod.size(); j++ ){
                                for(int i =0; i < r1prod.get(j).size(); i ++){      //loop to look through all the reactions involved production side of the rate function
                                    if(r1prod.get(j).get(i).equals(colname)){       //loop through all the PRODUCTS involved in the reaction in the previous loop                                        
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(j));
                                        
                                        if(isModInCol(mod2reactionmap, rowMRF.getProduction().get(j), colname)){
                                            if(isitDoubleModifier(rowMRF.getProduction().get(j))==false){
                                            prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivateWithModSubProd(false, i);
                                            prodreactions+=(int) stoichio1.get(j)+"Pmp"+(i+1);
                                            }else{
                                                prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivateWithDoubleModSubProd(isitActivatorOrInhibitor(rowMRF.getProduction().get(j), colname), false, i);
                                                prodreactions+=(int) stoichio1.get(j)+"PDmp"+(i+1);
                                            }
                                        }else{
                                            prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivate(false, i);
                                            prodreactions+=(int) stoichio1.get(j)+"Pp"+(i+1);
                                        }
                                    }
                                }
                            }
                        }
                        
                        if(inR2m==false&&inR2==false){    //if it is not involved in consumption reaction (production - 0)
                            consumereaction+=0;
                            consumereactions+=0;
                        }//now search whether its a sub/prod/mod/mod+sub/mod+prod in the consumption reaction 
                        
                        else if(inR2m==true&&inR2==false){              // just mod 
                            
                            for(int i =0;i<rowMRF.getConsumption().size();i++){
                                
                                if(isModInCol(mod2reactionmap, rowMRF.getConsumption().get(i), colname)){
                                    CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(i));
                                    if(isitDoubleModifier(rowMRF.getConsumption().get(i))==false){
                                    consumereaction+= (int) stoichio2.get(i) * function.getPartialDerivateWithMod();
                                    consumereactions+="justmod";
                                    }else{
                                        consumereaction+= (int) stoichio2.get(i) * function.getPartialDerivateWithDoubleMod(isitActivatorOrInhibitor(rowMRF.getConsumption().get(i), colname));
                                        consumereactions+="justDmod";
                                    }
                                }
                            }                            
                        }
                        
                        else if(inR2m==false&&inR2==true){             // just sub/prod
                            for(int j=0; j < r2sub.size(); j++ ){
                                for(int i =0; i < r2sub.get(j).size(); i ++){       //loop to look through all the reactions involved consumption side of the rate function
                                    if(r2sub.get(j).get(i).equals(colname)){        //loop through all the SUBSTRATES involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(j));
                                        consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivate(true, i);
                                        consumereactions+=(int) stoichio2.get(j) +"Cs"+(i+1);
                                    }
                                }
                            }
                            for(int j=0; j < r2prod.size(); j++ ){
                                for(int i =0; i < r2prod.get(j).size(); i ++){      //loop to look through all the reactions involved consumption side of the rate function
                                    if(r2prod.get(j).get(i).equals(colname)){       //loop through all the PRODUCTS involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(j));
                                        consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivate(false, i);
                                        consumereactions+=(int) stoichio2.get(j) +"Cp"+(i+1);
                                    }
                                }
                            }
                        }
                        
                        else if(inR2m==true&&inR2==true){              // mod +sub/prod
                            for(int j=0; j < r2sub.size(); j++ ){
                                for(int i =0; i < r2sub.get(j).size(); i ++){       //loop to look through all the reactions involved consumption side of the rate function
                                    if(r2sub.get(j).get(i).equals(colname)){        //loop through all the SUBSTRATES involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(j));
                                        if(isModInCol(mod2reactionmap, rowMRF.getConsumption().get(j), colname)){
                                            if(isitDoubleModifier(rowMRF.getConsumption().get(j))==false){
                                            consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivateWithModSubProd(true, i);
                                            consumereactions+=(int) stoichio2.get(j) +"Cms"+(i+1);
                                            }else{
                                                consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivateWithDoubleModSubProd(isitActivatorOrInhibitor(rowMRF.getConsumption().get(j), colname), true, i);
                                                consumereactions+=(int) stoichio2.get(j) +"CDms"+(i+1);
                                            }
                                        }else{
                                            consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivate(true, i);
                                            consumereactions+=(int) stoichio2.get(j) +"Cs"+(i+1);
                                        }
                                    }
                                }
                            }
                            for(int j=0; j < r2prod.size(); j++ ){
                                for(int i =0; i < r2prod.get(j).size(); i ++){      //loop to look through all the reactions involved consumption side of the rate function
                                    if(r2prod.get(j).get(i).equals(colname)){       //loop through all the PRODUCTS involved in the reaction in the previous loop
                                        CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(j));
                                        if(isModInCol(mod2reactionmap, rowMRF.getConsumption().get(j), colname)){
                                            if(isitDoubleModifier(rowMRF.getConsumption().get(j))==false){
                                            consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivateWithModSubProd(false, i);
                                            consumereactions+=(int) stoichio2.get(j) +"Cmp"+(i+1);
                                            }else{
                                                consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivateWithDoubleModSubProd(isitActivatorOrInhibitor(rowMRF.getConsumption().get(j), colname), false, i);
                                                consumereactions+=(int) stoichio2.get(j) +"CDmp"+(i+1);
                                            }
                                        }else{
                                            consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivate(false, i);
                                            consumereactions+=(int) stoichio2.get(j) +"Cp"+(i+1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                overall=prodreaction-consumereaction;
                overalls += prodreactions+"-"+consumereactions;
                matrix[row][col]=overall;
                stringmatrix[row][col]=overalls;
            }
        }
        
        Jmatrix = MatrixUtils.createRealMatrix(matrix);
        
    }
    
    private HashMap modifier2reaction(ArrayList<ModelReaction> Reactions){
        HashMap mod2reaction = new HashMap<>();
        
        for(ModelReaction reaction: Reactions){
            if(reaction.getRegulation()>1){
                mod2reaction.put(reaction.getName(), reaction.getModifier());
            }else{
                ArrayList NA = new ArrayList();
                NA.add(new Compound("N/A", 0.0));
                mod2reaction.put(reaction.getName(), NA);
            }            
        }
        
        return mod2reaction;
    }
    
    public RealMatrix getMatrix(){
        return Jmatrix;
    }
    
    public int getDimension(){
        return numComp;
    }
    
    public double[][] getM(){
        return matrix;
    }
    
    public String[][] getMm(){
        return stringmatrix;
    }
    
    private boolean isitDoubleModifier(String ReactionName){
        boolean isit = false;
        for(ModelReaction reactions : Reactions){
            if(reactions.getName().equals(ReactionName)){
                if(reactions.getRegulation()==4){
                    
                    isit=true;
                }
            }
        }
        return isit;
    }
    
    private int isitActivatorOrInhibitor(String ReactionName, String ColName){
        int activatorORinhibitor = 0; // 1 for activator 2 for inhibitor;
        for(ModelReaction reactions : Reactions){
            if(reactions.getName().equals(ReactionName)){
                String modname = reactions.getModifier().get(0).getName();
                String modname2 = reactions.getModifier().get(1).getName();
                if(modname.equals(ColName)){
                    activatorORinhibitor=1;
                }
                if(modname2.equals(ColName)){
                    activatorORinhibitor=2;
                }
            }
        }
        return activatorORinhibitor;
    }
    
    private boolean isModInCol(HashMap modreactionmap, String ReactionName, String ColName){
        boolean inCol = false;
        
        ArrayList<Compound> modifier = (ArrayList) modreactionmap.get(ReactionName);
        for(Compound comp : modifier){
            if(comp.getName().equals(ColName)){
                inCol=true;
            }
        }
                
        return inCol;
    }
}
