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
    private double[] uniUniParameters = new double[]{1, 1, 1, 1};
    private double[] uniBiRandomParameters = new double[]{1, 1, 1, 1, 1, 1};
    private double[] biUniRandomParameters = new double[]{1, 1, 1, 1, 1, 1};
    private double[] biBiRandomParameters = new double[]{1, 1, 1, 1, 1, 1, 1, 1};
    private boolean CKinetics = true, proteinkinetics; 
    private String[] subName, subID, prodName, prodID;
    private ArrayList<Double> CKparameters = new ArrayList<>();  //parameters for convenience kinetics
    private int regulation;                     //determines type of regulation, 1 for none, 2 for activation, 3 for inhibition
    private ArrayList substoichio, prodstoichio;
    private Compound modifier;
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
            modifier=mReaction.getModifier();
        }
        
        if (CKinetics==true){
            for (int i = 0; i < substrates.size(); i++) {
                subID[i] = substrates.get(i).getID();
                subName[i] = substrates.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }
            
            for (int i = 0; i < products.size(); i++) {
                prodID[i] = products.get(i).getID();
                prodName[i] = products.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }
        }else{
            a = substrates.get(0).getID();            
            aName = substrates.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            if (substrates.size() > 1) {
                b = substrates.get(1).getID();
                bName = substrates.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }

            p = products.get(0).getID();
            pName = products.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            if (products.size() > 1) {
                q = products.get(1).getID();
                qName = products.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }
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
            modifier=mReaction.getModifier();
        }
        
        if (CKinetics==true){
            for (int i = 0; i < substrates.size(); i++) {
                subID[i] = substrates.get(i).getID();
                subName[i] = substrates.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }
            
            for (int i = 0; i < products.size(); i++) {
                prodID[i] = products.get(i).getID();
                prodName[i] = products.get(i).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }
        }else{
            a = substrates.get(0).getID();            
            aName = substrates.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            if (substrates.size() > 1) {
                b = substrates.get(1).getID();
                bName = substrates.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }

            p = products.get(0).getID();
            pName = products.get(0).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            if (products.size() > 1) {
                q = products.get(1).getID();
                qName = products.get(1).getName().replaceAll("\\W", "").replaceAll("\\s", "");
            }
        }

        e = enzyme.getID().replaceAll("\\W", "").replaceAll("\\s", "");
    }
    
    public KineticLaw getKineticLaw() {

        String reactionRate = "";

        int nSub = substrates.size();
        int nProds = products.size();
        if (CKinetics == true){
            reactionRate = getCKinetics();
        }else{
            if (nSub == 1 && nProds == 1) {
                reactionRate = get11Reaction();

            } else if (nSub == 2 && nProds == 1) {
                reactionRate = get21Reaction();

            } else if (nSub == 1 && nProds == 2) {
                reactionRate = get12Reaction();

            } /*else if (nSub == 2 && nProds == 2) {
                reactionRate = get22Reaction();
            }*/
            else{
                reactionRate = get22Reaction();
            }
        }
        
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
    
    public double[] getUniUniParameters() {
        return uniUniParameters;
    }

    public double[] getUniBiRandomParameters() {
        return uniBiRandomParameters;
    }

    public double[] getBiUniRandomParameters() {
        return biUniRandomParameters;
    }

    public double[] getBiBiRandomParameters() {
        return biBiRandomParameters;
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
        
        if (regulation == 2){
            Ka.setId("Ka"+modifier.getName().replaceAll("\\W", "").replaceAll("\\s", ""));
            parameters.add(Ka);
        }else if(regulation == 3){
            Ki.setId("Ki"+modifier.getName().replaceAll("\\W", "").replaceAll("\\s", ""));
            parameters.add(Ki);
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
        if (regulation==1){
            top +=    "<apply>\n" +
                        "<times/>\n" +
                        "<ci>" + e+ "  </ci>\n" +
                        "<apply>" + "\n"+
                        "<minus/>\n"+
                        sub+prod+
                         "</apply>" + "\n"+
                        "</apply>" + "\n";
        
        } else if (regulation == 2){
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
                    "          <ci>"+ modifier.getID() +"</ci>\n" +
                    "          <ci> Ka"+modifier.getName().replaceAll("\\W", "").replaceAll("\\s", "") +"</ci>\n" + 
                    "        </apply>\n" +
                    "      </apply>\n" +
                    "    </apply>\n" +
                    "<apply>" + "\n"+
                    "<minus/>\n"+
                        sub+prod+
                         "</apply>" + "\n"+
                        "</apply>" + "\n";
        }else {
            top+=   "<apply>\n" +
                    "    <times/>\n" +
                    "    <apply>\n" +
                    "      <times/>\n" +
                    "<ci>" + e+ "  </ci>\n" +
                    "<apply>\n" +
                    "        <divide/>\n" +
                    "        <ci>  Ki"+modifier.getName().replaceAll("\\W", "").replaceAll("\\s", "") +" </ci>\n" +
                    "        <apply>\n" +
                    "          <plus/>\n" +
                    "          <ci>  Ki"+modifier.getName().replaceAll("\\W", "").replaceAll("\\s", "")+" </ci>\n" +     
                    "          <ci> "+ modifier.getID() +" </ci>\n"+      
                    "        </apply>\n" +
                    "      </apply>\n" +
                    "    </apply>\n" +
                    "<apply>" + "\n"+
                    "<minus/>\n"+
                    sub+prod+
                     "</apply>" + "\n"+
                    "</apply>" + "\n";
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
    
    private String get11Reaction() {
        Parameter KA = new Parameter("K" + aName);
        KA.setValue(getUniUniParameters()[0]);
        Parameter KB = new Parameter("K" + pName);
        KB.setValue(getUniUniParameters()[1]);
        Parameter KmA = new Parameter("Km" + aName);
        KmA.setValue(getUniUniParameters()[2]);
        Parameter KmP = new Parameter("Km" + pName);
        KmP.setValue(getUniUniParameters()[3]);

        parameters.add(KA);
        parameters.add(KB);
        parameters.add(KmA);
        parameters.add(KmP);
        return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" \n >"
                + "<apply>" + "\n"
                + "   <divide />" + "\n"
                + "   <apply>" + "\n"
                + "       <minus />" + "\n"
                + "       <apply>" + "\n"
                + "           <times />" + "\n"
                + "           <ci>K" + aName + "</ci>" + "\n"
                + "           <ci>" + e + "</ci>" + "\n"
                + "           <ci>" + a + "</ci>" + "\n"
                + "       </apply>" + "\n"
                + "       <apply>" + "\n"
                + "           <times />" + "\n"
                + "           <ci>K" + pName + "</ci>" + "\n"
                + "           <ci>" + e + "</ci>" + "\n"
                + "           <ci>" + p + "</ci>" + "\n"
                + "       </apply>" + "\n"
                + "   </apply>" + "\n"
                + "       <apply>" + "\n"
                + "           <plus />" + "\n"
                + "           <cn type=\"integer\">1</cn>" + "\n"
                + "       <apply>" + "\n"
                + "           <divide />" + "\n"
                + "           <ci>" + a + "</ci>" + "\n"
                + "           <ci>Km" + aName + "</ci>" + "\n"
                + "       </apply>" + "\n"
                + "       <apply>" + "\n"
                + "           <divide />" + "\n"
                + "           <ci>" + p + "</ci>" + "\n"
                + "           <ci>Km" + pName + "</ci>" + "\n"
                + "       </apply>" + "\n"
                + "   </apply>" + "\n"
                + "</apply>" + "\n"
                + "</math>";
    }
    
    private String get21Reaction() {
        
        Parameter KiA = new Parameter("Ki" + aName);
        KiA.setValue(getBiUniRandomParameters()[2]);
        Parameter KiB = new Parameter("Ki" + bName);
        KiB.setValue(getBiUniRandomParameters()[4]);
        Parameter KmB = new Parameter("Km" + bName);
        KmB.setValue(getBiUniRandomParameters()[5]);
        Parameter KmP = new Parameter("Km" + pName);
        KmP.setValue(getBiUniRandomParameters()[3]);
        Parameter Vr = new Parameter("Vr" + e);
        Vr.setValue(getBiUniRandomParameters()[1]);
        Parameter Vf = new Parameter("Vf" + e);
        Vf.setValue(getBiUniRandomParameters()[0]);

        parameters.add(Vf);
        parameters.add(Vr);
        parameters.add(KiA);
        parameters.add(KmB);
        parameters.add(KiB);
        parameters.add(KmP);
        
        return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" \n> "         
            + "<apply>" + "\n"
            + "  <divide/>" + "\n" 
            + " <apply>" + "\n"
            + "    <minus/>" + "\n" 
            + "    <apply>" + "\n"
            + "      <divide/>" + "\n"
            + "      <apply>" + "\n"
            + "        <times/>" + "\n"
            + "        <apply>" + "\n"
            + "          <times/>" + "\n"
            + "          <apply>" + "\n"
            + "            <times/>" + "\n"
            + "            <ci> Vf" + e+ "</ci>" + "\n"
            + "            <ci>" + e + "</ci>" + "\n"
            + "          </apply>" + "\n"
            + "          <ci>" + a + "</ci>" + "\n"
            + "        </apply>" + "\n"
            + "        <ci>" + b + "</ci>" + "\n"
            + "      </apply>" + "\n"
            + "      <apply>" + "\n"
            + "        <times/>" + "\n"
            + "        <ci>Ki" + aName + "</ci>" + "\n"
            + "        <ci>Km" + bName + "</ci>" + "\n"
            + "      </apply>" + "\n"
            + "    </apply>" + "\n"
            + "    <apply>" + "\n"
            + "      <divide/>" + "\n"
            + "      <apply>" + "\n"
            + "        <times/>" + "\n"
            + "        <apply>" + "\n"
            + "          <times/>" + "\n"
            + "          <ci> Vr" + e + "  </ci>" + "\n"
            + "          <ci> " + e + " </ci>" + "\n"
            + "        </apply>" + "\n"
            + "        <ci> " + p + " </ci>" + "\n"
            + "      </apply>" + "\n"
            + "      <ci> Km" + pName + " </ci>" + "\n"
            + "    </apply>" + "\n"
            + "  </apply>" + "\n"
            + "  <apply>" + "\n"
            + "    <plus/>" + "\n"
            + "    <apply>" + "\n"
            + "      <plus/>" + "\n"
            + "      <apply>" + "\n"
            + "        <plus/>" + "\n"
            + "        <apply>" + "\n"
            + "          <plus/>" + "\n"
            + "          <cn type=\"integer\"> 1 </cn>" + "\n"
            + "          <apply>" + "\n"
            + "            <divide/>" + "\n"
            + "               <ci>" + a + "</ci>" + "\n"
            + "               <ci>Ki" + aName + "</ci>" + "\n"
            + "          </apply>" + "\n"
            + "        </apply>" + "\n"
            + "        <apply>" + "\n"
            + "          <divide/>" + "\n"
            + "               <ci>" + b + "</ci>" + "\n"
            + "               <ci>Ki" + bName + "</ci>" + "\n"
            + "        </apply>" + "\n"
            + "      </apply>" + "\n"
            + "      <apply>" + "\n"
            + "        <divide/>" + "\n"
            + "        <apply>" + "\n"
            + "          <times/>" + "\n"
            + "               <ci>" + a + "</ci> " + "\n"
            + "               <ci>" + b + "</ci>" + "\n"
            + "        </apply>" + "\n"
            + "        <apply>" + "\n"
            + "          <times/>" + "\n"
            + "               <ci>Ki" + aName + "</ci> " + "\n"
            + "               <ci>Km" + bName + "</ci> " + "\n"
            + "        </apply>" + "\n"
            + "      </apply>" + "\n"
            + "    </apply>" + "\n"
            + "    <apply>" + "\n"
            + "      <divide/>" + "\n"
            + "               <ci>" + p + "</ci>" + "\n"
            + "               <ci>Km" + pName + "</ci>" + "\n"
            + "    </apply>" + "\n"
            + "  </apply>" + "\n"
            + "</apply>" + "\n"
            + "</math>" + "\n";
    }
    
    private String get12Reaction() {
        
        Parameter KmA = new Parameter("Km" + aName);
        KmA.setValue(getUniBiRandomParameters()[2]);
        Parameter KmP = new Parameter("Km" + pName);
        KmP.setValue(getUniBiRandomParameters()[4]);
        Parameter KiQ = new Parameter("Ki" + qName);
        KiQ.setValue(getUniBiRandomParameters()[5]);
        Parameter Vr = new Parameter("Vr" + e);
        Vr.setValue(getUniBiRandomParameters()[1]);
        Parameter KiP = new Parameter ("Ki" + pName);
        KiP.setValue(getUniBiRandomParameters()[3]);
        Parameter Vf = new Parameter("Vf" + e);
        Vf.setValue(getUniBiRandomParameters()[0]);
        
        parameters.add(Vf);
        parameters.add(Vr);
        parameters.add(KmA);
        parameters.add(KiP);
        parameters.add(KmP);
        parameters.add(KiQ);
        
        return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" \n >"
                + "<apply>" + "\n"
                + "<divide /> " + "\n"
                + "<apply>" + "\n"
                + "<minus /> " + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>Vf" + e + "</ci> " + "\n"
                + "<ci>" + e + "</ci> " + "\n"
                + "<ci>" + a + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "<ci>Km" + aName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>Vr" + e +"</ci>" + "\n"
                + "<ci>" + e + "</ci> " + "\n"
                + "<ci>" + p + "</ci>" + "\n"
                + "<ci>" + q + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>Km" + pName + "</ci>" + "\n"
                + "<ci>Ki" + qName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<plus /> " + "\n"
                + "<cn type=\"integer\">1</cn>" + "\n"
                + "<apply>" + "\n"
                + " <divide /> " + "\n"
                + "<ci>" + a + "</ci>" + "\n"
                + "<ci>Km" + aName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<ci>" + p + "</ci>" + "\n"
                + "<ci>Ki" + pName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + " <times /> " + "\n"
                + "<ci>" + p + "</ci>" + "\n"
                + "<ci>" + q + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>Km" + pName + "</ci> " + "\n"
                + "<ci>Ki" + qName + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide /> " + "\n"
                + "<ci>" + q + "</ci> " + "\n"
                + "<ci>Ki" + qName + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + " </math>" + "\n";
        
    }
    
    private String get22Reaction() {
        
        Parameter Vf = new Parameter("Vf" + e);
        Vf.setValue(getBiBiRandomParameters()[0]);
        Parameter Vr = new Parameter("Vr" + e);
        Vr.setValue(getBiBiRandomParameters()[1]);
        Parameter KiA = new Parameter("Ki" + aName);
        KiA.setValue(getBiBiRandomParameters()[2]);
        Parameter KmB = new Parameter("Km" + bName);
        KmB.setValue(getBiBiRandomParameters()[3]);
        Parameter KiB = new Parameter("Ki" + bName);
        KiB.setValue(getBiBiRandomParameters()[4]);
        Parameter KmP = new Parameter("Km" + pName);
        KmP.setValue(getBiBiRandomParameters()[5]);
        Parameter KiP = new Parameter("Ki" + pName);
        KiP.setValue(getBiBiRandomParameters()[6]);
        Parameter KiQ = new Parameter("Ki" + qName);
        KiQ.setValue(getBiBiRandomParameters()[7]);
        
        parameters.add(Vf);
        parameters.add(Vr);
        parameters.add(KiA);
        parameters.add(KmB);
        parameters.add(KiB);
        parameters.add(KmP);
        parameters.add(KiP);
        parameters.add(KiQ);
        
        return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" \n >"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + "<minus />" + "\n"
                + "<apply>" + "\n"
                + "<divide /> " + "\n"
                + "<apply>" + "\n"
                + "<times /> " + "\n"
                + "<ci>Vf" + e + "</ci> " + "\n"
                + "<ci>" + e + "</ci> " + "\n"
                + "<ci>" + a + "</ci> " + "\n"
                + "<ci>" + b + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>" + "Ki" + aName + "</ci> " + "\n"
                + "<ci>" + "Km" + bName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>Vr" + e + "</ci>" + "\n"
                + "<ci>" + e + "</ci> " + "\n"
                + "<ci>" + p + "</ci> " + "\n"
                + "<ci>" + q + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>" + "Km" + pName + "</ci> " + "\n"
                + "<ci>" + "Ki" + qName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<plus /> " + "\n"
                + "<cn type=\"integer\">1</cn>"
                + "<apply>" + "\n"
                + "<divide /> " + "\n"
                + "<ci>" + a + "</ci>" + "\n"
                + "<ci>" + "Ki" + aName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<ci>" + b + "</ci> " + "\n"
                + "<ci>" + "Ki" + bName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide /> " + "\n"
                + "<ci>" + p + "</ci>" + "\n"
                + "<ci>" + "Ki" + pName + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<ci>" + q + "</ci>" + "\n"
                + "<ci>" + "Ki" + qName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + " <apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + "<times /> " + "\n"
                + "<ci>" + a + "</ci> " + "\n"
                + "<ci>" + b + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>" + "Ki" + aName + "</ci> " + "\n"
                + "<ci>" + "Km" + bName + "</ci> " + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<divide />" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>" + p + "</ci>" + "\n"
                + "<ci>" + q + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "<apply>" + "\n"
                + "<times />" + "\n"
                + "<ci>" + "Km" + pName + "</ci>" + "\n"
                + "<ci>" + "Ki" + qName + "</ci>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "</apply>" + "\n"
                + "</math>";
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
    
    
    public void setBiBiRandomParameters(double[] d) {
        biBiRandomParameters = d;
    }
    
    public void setUniUniParameters(double[] d) {
        uniUniParameters = d;
    }
    
    public void setUniBiRandomParameters(double[] d) {
        uniBiRandomParameters = d;
    }
    
    public void setBiUniRandomParameters(double[] d) {
        biUniRandomParameters = d;
    }
    
    public void setCKineticsParameters(double[] d){
        
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
    
    public Compound getModifier(){
        return modifier;
    }
    public void setRegulation(int reg){
        regulation=reg;
    }
    public void setModifier(Compound comp){
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
