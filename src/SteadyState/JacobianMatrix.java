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
    
    public JacobianMatrix (VariableVector Vvector, ArrayList<Compound> Compounds, ArrayList<ModelReaction> Reactions, double[] Parameters){
        
        FunctionBuilder allfunctions = new FunctionBuilder(Vvector, Reactions, Parameters);
        HashMap functionmap = allfunctions.returnFunctions();   //build the MM Function for all the enzymatic reactions
        HashMap mod2reactionmap = modifier2reaction(Reactions);
        numComp = Compounds.size();
        
        ArrayList<String> modifiers=new ArrayList<>();
        for(ModelReaction reaction : Reactions){
            if(reaction.getRegulation()>1){
                modifiers.add(reaction.getModifier().getName());
            }
        }
        
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
                ArrayList<CKFunction> productionFunction = new ArrayList<>();
                ArrayList<CKFunction> consumptionFunction = new ArrayList<>();
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
                
                for (String name : modifiers){
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
                        inR1 = false;
                        inR2 = false;
                        
                        for(int i =0;i<rowMRF.getProduction().size();i++){
                            if(mod2reactionmap.get(rowMRF.getProduction().get(i)).equals(colname)){
                                inR1m=true;
                            }
                        }
                        
                        for(int i =0;i<rowMRF.getConsumption().size();i++){
                            if(mod2reactionmap.get(rowMRF.getConsumption().get(i)).equals(colname)){
                                inR2m=true;
                            }
                        }
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
                        
                        if(inR1m==false&&inR1==false){    //if it is not involved in production it becomes 0-consumption
                            prodreaction+=0;
                            prodreactions+=0;
                        }  //now search whether its a sub/prod/mod/mod+sub/mod+prod in the production reaction 
                        
                        else if(inR1m==true&&inR1==false){      // just mod
                            prodreactions+="justmod";
                            for(int i =0;i<rowMRF.getProduction().size();i++){
                                if(mod2reactionmap.get(rowMRF.getProduction().get(i)).equals(colname)){
                                    CKFunction function = (CKFunction) functionmap.get(rowMRF.getProduction().get(i));
                                    prodreaction+= (int) stoichio1.get(i) * function.getPartialDerivateWithMod();
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
                                        if(mod2reactionmap.get(rowMRF.getProduction().get(j)).equals(colname)){
                                            prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivateWithModSubProd(true, i);
                                            prodreactions+=(int) stoichio1.get(j)+"Pms"+(i+1);
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
                                        if(mod2reactionmap.get(rowMRF.getProduction().get(j)).equals(colname)){
                                            prodreaction+= (int) stoichio1.get(j) * function.getPartialDerivateWithModSubProd(false, i);
                                            prodreactions+=(int) stoichio1.get(j)+"Pmp"+(i+1);
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
                            consumereactions+="justmod";
                            for(int i =0;i<rowMRF.getConsumption().size();i++){
                                if(mod2reactionmap.get(rowMRF.getConsumption().get(i)).equals(colname)){
                                    CKFunction function = (CKFunction) functionmap.get(rowMRF.getConsumption().get(i));
                                    consumereaction+= (int) stoichio2.get(i) * function.getPartialDerivateWithMod();
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
                                        if(mod2reactionmap.get(rowMRF.getConsumption().get(j)).equals(colname)){
                                            consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivateWithModSubProd(true, i);
                                            consumereactions+=(int) stoichio2.get(j) +"Cms"+(i+1);
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
                                        if(mod2reactionmap.get(rowMRF.getConsumption().get(j)).equals(colname)){
                                            consumereaction+= (int) stoichio2.get(j) * function.getPartialDerivateWithModSubProd(false, i);
                                            consumereactions+=(int) stoichio2.get(j) +"Cmp"+(i+1);
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
                mod2reaction.put(reaction.getName(), reaction.getModifier().getName());
            }else{
                mod2reaction.put(reaction.getName(), "N/A");
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
}