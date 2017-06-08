/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import Reaction.MetabolicReaction;
import Species.*;
import java.util.ArrayList;
import org.sbml.jsbml.ASTNode;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.LocalParameter;
import org.sbml.jsbml.Parameter;

/**
 *
 * @author chuanfuyap
 */
public class ModelReaction {
    
    private MetabolicReaction mReaction;
    private ArrayList<Compound> substrates;
    private ArrayList<Compound> products;
    private ArrayList<Parameter> parameters=new ArrayList();
    private ArrayList<LocalParameter> localparameters;
    private Enzyme enzyme;
    private String a, b, e, p, q;                       //for storing species id.
    private String aName, bName, pName, qName;          //for storing species Name.
    private boolean proteinkinetics;
    private String[] subName, subID, prodName, prodID;
    private ArrayList<Double> CKparameters = new ArrayList<>();  //parameters for convenience kinetics
    private int regulation;                     //determines type of regulation, 1 for none, 2 for activation, 3 for inhibition
    private ArrayList substoichio, prodstoichio;
    private ArrayList<Compound> modifier;
    private ArrayList<Parameter> PKparameters = new ArrayList();
    private ArrayList<Parameter> GKparameters = new ArrayList();
    
    public ModelReaction(Enzyme enzyme, MetabolicReaction mReaction) {
        
        this.enzyme = enzyme;
        this.mReaction = mReaction;
        substrates = mReaction.getSubstrates();
        products = mReaction.getProducts();
        parameters = new ArrayList<>();
        subName = new String[substrates.size()];
        subID = new String[substrates.size()];
        prodName = new String[products.size()];
        prodID = new String[products.size()];
        substoichio = mReaction.getSubstrateStoichio();
        prodstoichio = mReaction.getProductStoichio();
        proteinkinetics=mReaction.getPKinetics();
        regulation=mReaction.getRegulation();
        if (regulation>1){
            modifier=mReaction.getModifiers();
        }
        
        for (int i = 0; i < substrates.size(); i++) {
            subID[i] = substrates.get(i).getID();
            subName[i] = substrates.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
        }

        for (int i = 0; i < products.size(); i++) {
            prodID[i] = products.get(i).getID();
            prodName[i] = products.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
        }
        
        e = enzyme.getID().replaceAll("\\W", "").replaceAll("\\s", "");

    }
    
    public ModelReaction(Enzyme enzyme, MetabolicReaction mReaction, ArrayList<LocalParameter> Parameters) {

        this.enzyme = enzyme;
        this.mReaction = mReaction;
        substrates = mReaction.getSubstrates();
        products = mReaction.getProducts();
        localparameters = Parameters;
        subName = new String[substrates.size()];
        subID = new String[substrates.size()];
        prodName = new String[products.size()];
        prodID = new String[products.size()];
        substoichio = mReaction.getSubstrateStoichio();
        prodstoichio = mReaction.getProductStoichio();
        proteinkinetics=mReaction.getPKinetics();
        regulation=mReaction.getRegulation();
        if (regulation>1){
            modifier=mReaction.getModifiers();
        }
        
        for (int i = 0; i < substrates.size(); i++) {
            subID[i] = substrates.get(i).getID();
            subName[i] = substrates.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
        }

        for (int i = 0; i < products.size(); i++) {
            prodID[i] = products.get(i).getID();
            prodName[i] = products.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
        }
        
        e = enzyme.getID().replaceAll("\\W", "").replaceAll("\\s", "");
    }
    
    public KineticLaw getKineticLaw() {

        String reactionRate = "";

        int nSub = substrates.size();
        int nProds = products.size();
        reactionRate = getCKinetics();

        KineticLaw kl = new KineticLaw();
        ASTNode mathmlnode = ASTNode.readMathMLFromString(reactionRate);
        kl.setMath(mathmlnode);
        return kl;
    }
    
    public boolean getPkinetics(){
        return proteinkinetics;
    }
    
    public void setReactionID(String ID) {
        mReaction.setReactionId(ID);
    }

    public ArrayList<Parameter> getParameters() {
        return parameters;
    }
    
    public ArrayList<LocalParameter> getlocalparameters(){
        return localparameters;
    }
    
    public String getName() {
        return mReaction.getName();
    }
    
    public String getID() {
        return mReaction.getReactionID();
    }
        
    private String getCKinetics(){
        parameters.clear();
        String equation="";
        
        Parameter Vf = new Parameter("Vf" + e.replaceAll("\\W", "").replaceAll("\\s", ""));
        Vf.setValue(1);
        Parameter Vr = new Parameter("Vr" + e.replaceAll("\\W", "").replaceAll("\\s", ""));
        Vr.setValue(1);
        Parameter Ka = new Parameter();
        Ka.setValue(1);
        Parameter Ki = new Parameter();
        Ki.setValue(1);
        ArrayList<Parameter> subparams = new ArrayList<>();
        ArrayList<Parameter> prodparams = new ArrayList<>();
        parameters.add(Vf);
        parameters.add(Vr);
        for (int i = 0; i < substrates.size(); i++) {
            Parameter subparam = new Parameter("Km" + subName[i]);
            subparam.setValue(1);
            subparams.add(subparam);
        }
        for (int i = 0; i < products.size(); i++) {
            Parameter prodparam = new Parameter("Km" + prodName[i]);
            prodparam.setValue(1);
            prodparams.add(prodparam);
        }
        
        switch (regulation) {
            case 2:
                Ka.setId("Ka"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", ""));
                parameters.add(Ka);
                break;
            case 3:
                Ki.setId("Ki"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", ""));
                parameters.add(Ki);
                break;
            case 4:
                Ka.setId("Ka"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", ""));
                parameters.add(Ka);
                Ki.setId("Ki"+modifier.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", ""));
                parameters.add(Ki);
                break;
            default:
                break;
        }
        
        for (int i = 0; i < subparams.size(); i++) {
            parameters.add(subparams.get(i));
        }
        for (int i = 0; i < prodparams.size(); i++) {
            parameters.add(prodparams.get(i));
        }
        
        
        //building equation for the numerator
        String sub ="<apply>\n" +
                        "<times/>\n" +
                    "<ci>Vf" + e+ "  </ci>\n";
        if(substrates.size()==1){
            
            sub+=  "<apply>\n" +
                    "<power/>\n" +
                    "<apply>" + "\n"+
                    "<divide/>\n" +
                        "    <ci>"+ subID[0] +"</ci>\n" +
                        "    <ci> Km"+ subName[0] +"</ci>\n"+
                    "</apply>\n" +
                                "        <cn type=\"integer\">"+ (int)substoichio.get(0) +"</cn>\n" +
                    "</apply>" + "\n";
        }else{
            for(int i =0; i<substrates.size()-1;i++){
                sub += "<apply>" + "\n"
                    + "<times/>" + "\n";
            }
            for ( int i =0; i<substrates.size(); i++){
                
                sub+=    "<apply>\n" +
                                "        <power/>\n" +
                        "<apply>" + "\n"+
                        "<divide/>\n" +
                        "    <ci>"+ subID[i] +"</ci>\n" +
                        "    <ci> Km"+ subName[i] +"</ci>\n"+
                        "        </apply>\n" +
                                "        <cn type=\"integer\">"+ (int)substoichio.get(i) +"</cn>\n" +
                        "</apply>" + "\n";
            }
            for(int i =0; i<substrates.size()-1;i++){
                sub +="</apply>" + "\n";
            }
        }
        sub +="</apply>" + "\n";
        
        String prod="<apply>\n" +
                        "<times/>\n" +
                    "<ci>Vr" + e+ "  </ci>\n";
        if(products.size()==1){
            
            prod+=   "<apply>\n" +
                                "        <power/>\n"
                                        +"<apply>" + "\n"+
                        "<divide/>\n" +
                        "    <ci>"+ prodID[0] +"</ci>\n" +
                        "    <ci> Km"+ prodName[0] +"</ci>\n"+
                                        "        </apply>\n" +
                                "        <cn type=\"integer\">"+ (int) prodstoichio.get(0) +"</cn>\n" +
                    "</apply>" + "\n";
        }else{
            for(int i =0; i<products.size()-1;i++){
                prod += "<apply>" + "\n"
                    + "<times/>" + "\n";
            }
            
            for ( int i =0; i<products.size(); i++){
                prod+=   "<apply>\n" +
                                "        <power/>\n"+
                        "<apply>" + "\n"+
                        "<divide/>\n" +
                        "    <ci>"+ prodID[i] +"</ci>\n" +
                        "    <ci> Km"+ prodName[i] +"</ci>\n"+
                        "        </apply>\n" +
                                "        <cn type=\"integer\">"+ (int) prodstoichio.get(i)  +"</cn>\n" +
                        "</apply>" + "\n";
            }
            for(int i =0; i<products.size()-1;i++){
                prod+="</apply>" + "\n";
            }
        }
        prod+="</apply>" + "\n";
        String top ="";
        switch (regulation) {
            case 1:
                top +=  "<apply>\n" +
                        "<times/>\n" +
                        "<ci>" + e+ "  </ci>\n" +
                        "<apply>" + "\n"+
                        "<minus/>\n"+
                        sub+prod+
                        "</apply>" + "\n"+
                        "</apply>" + "\n";
                break;
            case 2:
                top+=   "<apply>\n" +
                        "    <times/>\n" +
                        "    <apply>\n" +
                        "      <times/>\n" +
                        "<ci>" + e+ "  </ci>\n" +
                        "      <apply>\n" +
                        "        <plus/>\n" +
                        "        <cn type=\"integer\"> 1 </cn>\n" +
                        "        <apply>\n" +
                        "          <divide/>\n" +
                        "          <ci>"+ modifier.get(0).getID() +"</ci>\n" +
                        "          <ci> Ka"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "") +"</ci>\n" +
                        "        </apply>\n" +
                        "      </apply>\n" +
                        "    </apply>\n" +
                        "<apply>" + "\n"+
                        "<minus/>\n"+
                        sub+prod+
                        "</apply>" + "\n"+
                        "</apply>" + "\n";
                break;
            case 3:
                top+=   "<apply>\n" +
                        "    <times/>\n" +
                        "    <apply>\n" +
                        "      <times/>\n" +
                        "<ci>" + e+ "  </ci>\n" +
                        "<apply>\n" +
                        "        <divide/>\n" +
                        "        <ci>  Ki"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "") +" </ci>\n" +
                        "        <apply>\n" +
                        "          <plus/>\n" +
                        "          <ci>  Ki"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "")+" </ci>\n" +
                        "          <ci> "+ modifier.get(0).getID() +" </ci>\n"+
                        "        </apply>\n" +
                        "      </apply>\n" +
                        "    </apply>\n" +
                        "<apply>" + "\n"+
                        "<minus/>\n"+
                        sub+prod+
                        "</apply>" + "\n"+
                        "</apply>" + "\n";
                break;
            case 4:
                top+=   "              <apply>\n" +
                        "                <times/>\n" +
                        "                <apply>\n" +
                        "                  <times/>\n" +
                        "                  <apply>\n" +
                        "                    <times/>\n" +
                        "                    <ci>" + e+ "  </ci>\n" +
                        "                    <apply>\n" +
                        "                      <plus/>\n" +
                        "                      <cn type=\"integer\"> 1 </cn>\n" +
                        "                      <apply>\n" +
                        "                        <divide/>\n" +
                        "                        <ci>"+ modifier.get(0).getID() +"</ci>\n" +
                        "                        <ci> Ka"+modifier.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "") +"</ci>\n" +
                        "                      </apply>\n" +
                        "                    </apply>\n" +
                        "                  </apply>\n" +
                        "                  <apply>\n" +
                        "                    <divide/>\n" +
                        "                    <ci>  Ki"+modifier.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", "") +" </ci>\n" +
                        "                    <apply>\n" +
                        "                      <plus/>\n" +
                        "                      <ci>  Ki"+modifier.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", "")+" </ci>\n" +
                        "                      <ci> "+ modifier.get(1).getID() +" </ci>\n"+
                        "                    </apply>\n" +
                        "                  </apply>\n" +
                        "                </apply>\n" +
                        "                <apply>\n" +
                        "                  <minus/>\n" +
                        sub+prod+
                        "</apply>" + "\n"+
                        "</apply>" + "\n";
                break;
            default:
                break;
        }
        
        //building equation for the denominator
        String btmsub ="";
        if(substrates.size()==1){
            if((int)substoichio.get(0)==1){
                btmsub+="  <apply>\n" +
                        "    <plus/>\n" +
                        "    <cn type=\"integer\"> 1 </cn>\n" +
                        "    <apply>\n" +
                        "      <divide/>\n" +
                        "    <ci>"+ subID[0] +"</ci>\n" +
                        "    <ci> Km"+ subName[0] +"</ci>\n"+
                        "    </apply>\n" +
                        "  </apply>\n";
            }else{
                for(int j=0; j<(int)substoichio.get(0);j++){
                            btmsub+=    "  <apply>\n" +
                                        "    <plus/>\n";
                }

                btmsub+=    "    <cn type=\"integer\"> 1 </cn>\n"+
                            "       <apply>\n" +
                            "          <divide/>\n" +
                            "          <ci>"+ subID[0] +"</ci>\n" +
                            "          <ci> Km"+ subName[0] +"</ci>\n"+
                            "        </apply>\n"+
                            "</apply>\n";
                for(int j=1; j<(int)substoichio.get(0);j++){
                    int count = 1;
                    count +=j;
                    btmsub+=    "<apply>\n" +
                                "        <power/>\n" +
                                "        <apply>\n" +
                                "          <divide/>\n" +
                                "    <ci>"+ subID[0] +"</ci>\n" +
                                "    <ci> Km"+ subName[0] +"</ci>\n"+
                                "        </apply>\n" +
                                "        <cn type=\"integer\">"+ count +"</cn>\n" +
                                "      </apply>\n"+
                                "</apply>\n";
                }
            }
        }else{
            for(int i =0; i<substrates.size()-1;i++){
                btmsub += "<apply>" + "\n"
                    + "<times/>" + "\n";
            }
            for ( int i =0; i<substrates.size(); i++){
                if((int)substoichio.get(i)==1){
                btmsub+="  <apply>\n" +
                        "    <plus/>\n" +
                        "    <cn type=\"integer\"> 1 </cn>\n" +
                        "    <apply>\n" +
                        "      <divide/>\n" +
                        "    <ci>"+ subID[i] +"</ci>\n" +
                        "    <ci> Km"+ subName[i] +"</ci>\n"+
                        "    </apply>\n"+
                        "   </apply>\n";
                }else{

                    for(int j=0; j<(int)substoichio.get(i);j++){
                                btmsub+=    "  <apply>\n" +
                                            "    <plus/>\n";
                    }

                    btmsub+=    "    <cn type=\"integer\"> 1 </cn>\n"+
                                "       <apply>\n" +
                                "          <divide/>\n" +
                                "          <ci>"+ subID[i] +"</ci>\n" +
                                "          <ci> Km"+ subName[i] +"</ci>\n"+
                                "        </apply>\n"+
                                "</apply>\n";
                    
                    for(int j=1; j<(int)substoichio.get(i);j++){
                        int count = 1;
                        count +=j;
                        btmsub+=    "<apply>\n" +
                                    "        <power/>\n" +
                                    "        <apply>\n" +
                                    "          <divide/>\n" +
                                    "    <ci>"+ subID[i] +"</ci>\n" +
                                    "    <ci> Km"+ subName[i] +"</ci>\n"+
                                    "        </apply>\n" +
                                    "        <cn type=\"integer\">"+ count +"</cn>\n" +
                                    "      </apply>\n"+
                                    "</apply>\n";
                    }
                }
                
            }
            for(int i =0; i<substrates.size()-1;i++){
                btmsub +="</apply>" + "\n";
            }
        }
        
        String btmprod="";
        if(products.size()==1){
            if((int)prodstoichio.get(0)==1){
                btmprod+="  <apply>\n" +
                        "    <plus/>\n" +
                        "    <cn type=\"integer\"> 1 </cn>\n" +
                        "    <apply>\n" +
                        "      <divide/>\n" +
                        "    <ci>"+ prodID[0] +"</ci>\n" +
                        "    <ci> Km"+ prodName[0] +"</ci>\n"+
                        "    </apply>\n" +
                        "  </apply>\n";
            }else{
                
                for(int j=0; j<(int)prodstoichio.get(0);j++){
                            btmprod+=    "  <apply>\n" +
                                        "    <plus/>\n";
                }

                btmprod+=    "    <cn type=\"integer\"> 1 </cn>\n"+
                            "       <apply>\n" +
                            "          <divide/>\n" +
                            "          <ci>"+ prodID[0] +"</ci>\n" +
                            "          <ci> Km"+ prodName[0] +"</ci>\n"+
                            "        </apply>\n"+
                        "</apply>\n";
                
                for(int j=1; j<(int)prodstoichio.get(0);j++){
                    int count = 1;
                    count +=j;
                    btmprod+=    "<apply>\n" +
                                "        <power/>\n" +
                                "        <apply>\n" +
                                "          <divide/>\n" +
                                "    <ci>"+ prodID[0] +"</ci>\n" +
                                "    <ci> Km"+ prodName[0] +"</ci>\n"+
                                "        </apply>\n" +
                                "        <cn type=\"integer\">"+ count +"</cn>\n" +
                                "      </apply>\n"+
                                "</apply>\n";
                }
                
            }
        }else{
            for(int i =0; i<products.size()-1;i++){
                btmprod += "<apply>" + "\n"
                    + "<times/>" + "\n";
            }
            for ( int i =0; i<products.size(); i++){
                if((int)prodstoichio.get(i)==1){
                btmprod+="  <apply>\n" +
                        "    <plus/>\n" +
                        "    <cn type=\"integer\"> 1 </cn>\n" +
                        "    <apply>\n" +
                        "      <divide/>\n" +
                        "    <ci>"+ prodID[i] +"</ci>\n" +
                        "    <ci> Km"+ prodName[i] +"</ci>\n"+
                        "    </apply>\n" +
                        "  </apply>\n";
            }else{
                
                for(int j=0; j<(int)prodstoichio.get(i);j++){
                            btmprod+=    "  <apply>\n" +
                                        "    <plus/>\n";
                }

                btmprod+=    "    <cn type=\"integer\"> 1 </cn>\n"+
                            "       <apply>\n" +
                            "          <divide/>\n" +
                            "          <ci>"+ prodID[i] +"</ci>\n" +
                            "          <ci> Km"+ prodName[i] +"</ci>\n"+
                            "        </apply>\n"+
                            "</apply>\n";
                
                for(int j=1; j<(int)prodstoichio.get(i);j++){
                    int count = 1;
                    count +=j;
                    btmprod+=    "<apply>\n" +
                                "        <power/>\n" +
                                "        <apply>\n" +
                                "          <divide/>\n" +
                                "    <ci>"+ prodID[i] +"</ci>\n" +
                                "    <ci> Km"+ prodName[i] +"</ci>\n"+
                                "        </apply>\n" +
                                "        <cn type=\"integer\">"+ count +"</cn>\n" +
                                "      </apply>\n"+
                                "</apply>\n";
                }
                
            }
                
            }
            for(int i =0; i<products.size()-1;i++){
                btmprod +="</apply>" + "\n";
            }
        }
        
        String btm = "<apply>\n" +
                        "    <minus/>\n" +
                        "    <apply>\n" +
                        "      <plus/>\n"+
                        btmsub+btmprod+
                        "</apply>" + "\n"+
                        "<cn type=\"integer\">"+ 1 +"</cn>\n" +
                        "</apply>" + "\n";
                
        
        
        equation += "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n" +
                    "<apply>\n" +
                    "<divide/>\n" +
                    top+
                    btm+
                    "</apply>\n"+
                    "</math>";
        
        return equation;
    }
    
    public MetabolicReaction getMyReaction() {
        return mReaction;
    }
    
    public String getReactionID() {
        return mReaction.getReactionID();
    }
    
    public Enzyme getEnzyme() {
        return enzyme;
    }
    
    public ArrayList<Compound> getSubstrates() {
        return mReaction.getSubstrates();
    }

    public ArrayList<Compound> getProducts() {
        return mReaction.getProducts();
    }

    public ArrayList getSubstratesStoichio() {
        return mReaction.getSubstrateStoichio();
    }

    public ArrayList getProductsStoichio() {
        return mReaction.getProductStoichio();
    }
        
    public KineticLaw getProteinKinetics(){
        Parameter translate = new Parameter("Ktl" + e);
        translate.setValue(1);
        Parameter degradation = new Parameter("Kdeg" + e);
        degradation.setValue(1);
        
        PKparameters.add(translate);
        PKparameters.add(degradation);
        
        String kinetics =  " <math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n" +
                "  <apply>\n" +
                "    <minus/>\n" +
                "    <apply>\n" +
                "      <times/>\n" +
                "      <ci> Ktl" + e+" </ci>\n" +
                "      <ci> mRNA"+e+" </ci>\n" +
                "    </apply>\n" +
                "    <apply>\n" +
                "      <times/>\n" +
                "      <ci> Kdeg"+e+" </ci>\n" +
                "      <ci> "+e+" </ci>\n" +
                "    </apply>\n" +
                "  </apply>\n" +
                "</math>";
        
        KineticLaw kl = new KineticLaw();
        ASTNode mathmlnode = ASTNode.readMathMLFromString(kinetics);
        kl.setMath(mathmlnode);
        return kl;
    }
    
    public ArrayList<Parameter> getPKparameters() {
        return PKparameters;
    }
    
    public KineticLaw getGeneKinetics(){
        Parameter transcription = new Parameter("Ktc" + e);
        transcription.setValue(1);
        Parameter degradation = new Parameter("KdegRNA" + e);
        degradation.setValue(1);
        
        GKparameters.add(transcription);
        GKparameters.add(degradation);
        
        String kinetics =  " <math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n" +
                "  <apply>\n" +
                "    <minus/>\n" +
                "      <ci> Ktc" + e+" </ci>\n" +
                "    <apply>\n" +
                "      <times/>\n" +
                "      <ci> KdegRNA"+e+" </ci>\n" +
                "      <ci> mRNA"+e+" </ci>\n" +
                "    </apply>\n" +
                "  </apply>\n" +
                "</math>";
        
        KineticLaw kl = new KineticLaw();
        ASTNode mathmlnode = ASTNode.readMathMLFromString(kinetics);
        kl.setMath(mathmlnode);
        return kl;
    }
    
    public ArrayList<Parameter> getGKparameters() {
        return GKparameters;
    }
    
    public int getRegulation(){
        return regulation;
    }
    
    public ArrayList<Compound> getModifier(){
        return modifier;
    }
    public void setRegulation(int reg){
        regulation=reg;
    }
    public void setModifier(ArrayList<Compound> comp){
        modifier=comp;
    }
    
    public void setPKinetics(boolean yesno){
        proteinkinetics=yesno;
    }
    public void setEnzyme(Enzyme enz){
        enzyme=enz;
    }
    
    public void EditSubstratesProducts(Enzyme newenzyme, ArrayList<Compound> newsubstrates, ArrayList newstoichio1, ArrayList<Compound> newproducts,ArrayList newstoichio2 ){
        this.enzyme = newenzyme;
        substrates = newsubstrates;
        products = newproducts;
        subName = new String[substrates.size()];
        subID = new String[substrates.size()];
        prodName = new String[products.size()];
        prodID = new String[products.size()];
        substoichio = newstoichio1;
        prodstoichio = newstoichio2;
        parameters = new ArrayList<>();
        
        for (int i = 0; i < substrates.size(); i++) {
            subID[i] = substrates.get(i).getID();
            subName[i] = substrates.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
        }

        for (int i = 0; i < products.size(); i++) {
            prodID[i] = products.get(i).getID();
            prodName[i] = products.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
        }

        e = enzyme.getID().replaceAll("\\W", "").replaceAll("\\s", "");
        
    }
    
    public String toString() {

        String output = " ";
        String temp = "";
        String temp2 = "";
        String ID = mReaction.getReactionID();

        for (int i = 0; i < mReaction.getSubstrates().size(); i++) {
            
            if (mReaction.getSubstrates().size() == 1) {
                temp +=  "";
            }
            

            if (i == mReaction.getSubstrates().size() - 1) {
                temp +="(" + mReaction.getSubstrateStoicCo(i)+ ") " + mReaction.getSubstrates().get(i).getName();
            } else {
                temp += "(" + mReaction.getSubstrateStoicCo(i) + ") " + mReaction.getSubstrates().get(i).getName() + " + ";
            }
        }

        for (int i = 0; i < products.size(); i++) {

            if (i == mReaction.getProducts().size() - 1) {
                temp2 += "(" +mReaction.getProductStoicCo(i) + ") " + mReaction.getProducts().get(i).getName() + "";
            } else {
                temp2 +="(" + mReaction.getProductStoicCo(i) + ") " + mReaction.getProducts().get(i).getName() + "  + ";
            }

        }
        
        output = ID+": "+ temp + " <-> "+ temp2;
        return output;
    }
    
    public ArrayList getSubStoichio(){
        return substoichio; 
    }
    
    public ArrayList getProdStoichio(){
        return prodstoichio;
    }
}
