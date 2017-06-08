/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import Species.*;
import SBML.ImportSBML;
import Reaction.MetabolicReaction;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import javax.swing.JOptionPane;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.ListOf;
import org.sbml.jsbml.LocalParameter;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.ModifierSpeciesReference;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLException;
import org.sbml.jsbml.SBMLReader;
import org.sbml.jsbml.SBMLWriter;
import org.sbml.jsbml.Species;
import org.sbml.jsbml.SpeciesReference;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class SBMLFile {
    
    private final SBMLDocument sbmldoc;
    private final Model model;     
    private final ListOf<Reaction> reactions;
    private KineticLaw KL=null;
    private ListOf<Reaction> modelReaction;
    private ListOf<Species> listofspecies;
    private ArrayList<ModelReaction>  allthereactions = new ArrayList<>();
    private ArrayList<Compound> allthespecies = new ArrayList<Compound>();
    private ArrayList<Enzyme> alltheenzymes = new ArrayList<Enzyme>();
    private String modelName;
    private ImportSBML sbmlmodel;
    private ListOf<LocalParameter> LLP;
    private ArrayList<String> parameternames = new ArrayList<>();
    private ArrayList<String> reactionIDlist = new ArrayList<>();
    
    public SBMLFile(File file) throws XMLStreamException, IOException, ModelOverdeterminedException {
        
        sbmldoc = new SBMLDocument (SBMLReader.read(file));
        model = sbmldoc.getModel();
        sbmlmodel = new ImportSBML(file); 
        modelName = sbmlmodel.getModelName();
        modelReaction = sbmlmodel.getReactions();
        listofspecies = sbmlmodel.getSpecies();
        addReactionsToReactionList(modelReaction);
        reactions = (ListOf) model.getListOfReactions();
        addSpeciesToSpeciesList(listofspecies);
        
    }
    private void addSpeciesToSpeciesList(ListOf<Species> listofspecies){
        for (int i = 0; i < listofspecies.size(); i++) {
            String name = listofspecies.get(i).getName();
            String ID = listofspecies.get(i).getId();
            int SBO = listofspecies.get(i).getSBOTerm();
            Compound newCompound = new Compound(listofspecies.get(i));            
            String sub = ID.substring(0,1);
            String sub2="";
            if(ID.length()>1){
                sub2 = ID.substring(1,2);
            }
            if(SBO==252 || (sub.equals("E") && sub2.matches("\\d"))){
                
                alltheenzymes.add(new Enzyme(listofspecies.get(i)));
            }else if (SBO ==250){
            }else{
                
                allthespecies.add(newCompound);
            }
        }
        Object[][] speciesinfo = new Object[allthespecies.size()][4];
        Object[][] enzymeinfo = new Object[alltheenzymes.size()][3];
        
        for(int i =0; i<allthespecies.size(); i++){
            speciesinfo[i][0]=i+1;
            speciesinfo[i][1]=allthespecies.get(i).getName();
            speciesinfo[i][2]=allthespecies.get(i).getConcentration();
            speciesinfo[i][3]=allthespecies.get(i).getBoundaryCondition();
        }
        
        for(int i =0; i<alltheenzymes.size(); i++){
            enzymeinfo[i][0]=i+1;
            enzymeinfo[i][1]=alltheenzymes.get(i).getName();
            enzymeinfo[i][2]=alltheenzymes.get(i).getConcentration();
        }
        
    }
    private void addReactionsToReactionList(ListOf<Reaction> listofreactions){
        
        for (int i = 0; i < listofreactions.size(); i++){           //reads the model's list of reactions.
            Reaction reaction = listofreactions.get(i);             //extracts each reaction found in the list, go through them 1 by 1 in the For loop.
            MetabolicReaction mReaction = new MetabolicReaction(reaction);
            ListOf modifiers = (ListOf) reaction.getListOfModifiers();
            ListOf reactionsubstrates = (ListOf) reaction.getListOfReactants();
            ListOf reactionproducts = (ListOf) reaction.getListOfProducts();
            Enzyme e = null;
            ArrayList<LocalParameter> parameters = new ArrayList<LocalParameter>();
            int regulation =1;
            ArrayList<Compound> mod = new ArrayList();
            String isitPK = reaction.getName();
            String[] tokens = isitPK.split("\\s");
            
            if((tokens.length>1&&"transcription".equals(tokens[1]))||(tokens.length>1&&"translation".equals(tokens[1]))){
                String pkreaction = tokens[0];
                for(int j =0; j<allthereactions.size();j++){
                    if(allthereactions.get(j).getName().equals(pkreaction)){
                        allthereactions.get(j).setPKinetics(true);
                    }
                }
            }else{
                for (int j = 0; j <reactionsubstrates.size(); j++){
                    SpeciesReference sf = (SpeciesReference) reactionsubstrates.get(j);
                    Species substratespecies = sf.getSpeciesInstance();
                    Compound c = new Compound(substratespecies);
                    mReaction.addSubstrate(c, (int) sf.getStoichiometry());
                }

                for (int k = 0; k <reactionproducts.size(); k++){
                    SpeciesReference sf = (SpeciesReference) reactionproducts.get(k);
                    Species productspecies = sf.getSpeciesInstance();
                    Compound c = new Compound(productspecies);
                    mReaction.addProduct(c, (int) sf.getStoichiometry());
                }

                KineticLaw KL=reaction.getKineticLaw();
                LLP=KL.getListOfLocalParameters();

                for (int count = 0; count < LLP.size(); count++){
                    parameters.add(LLP.get(count));
                }

                if (modifiers.size()==0){
                    Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
                    enzymeSpecies.setName("EnZ"+reaction.getName());
                    enzymeSpecies.setInitialConcentration(1);
                    e = new Enzyme(enzymeSpecies);
                    addEnzymesFromReaction(e);
                }else{
                    //sbo 20/206/207/597 are inhibitors
                    //sbo 461/535/534/533 are activators;
                    //sbo 460 is enzyme catalyst
                    boolean old = true;
                    for ( int l=0;l<modifiers.size();l++){
                        int sbo = modifiers.get(l).getSBOTerm();

                        if(sbo!=-1){
                            old=false;
                        }
                    }
                    if(old==true){
                        e = oldSBML( mReaction,  reaction,  modifiers, e);
                    }else{
                        for ( int l=0;l<modifiers.size();l++){
                            int sbo = modifiers.get(l).getSBOTerm();
                            if(sbo==460){
                                ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
                                Species s = msr.getSpeciesInstance();
                                e = new Enzyme(s);
                            }else if (sbo==461||sbo==535||sbo==534||sbo==533){
                                regulation = 2; 
                                ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
                                Species s = msr.getSpeciesInstance();
                                mod.add( new Compound(s) );
                            }else if (sbo==20||sbo==206||sbo==207||sbo==597){
                                regulation = 3; 
                                ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
                                Species s = msr.getSpeciesInstance();
                                mod.add( new Compound(s) );
                            }
                        }
                    }
                }

                if (regulation>1){
                    mReaction.setEnzyme(e);
                    mReaction.setModifier(mod);
                    if(mod.size()==2){
                        mReaction.setRegulation(4);
                    }else if (mod.size()==1){
                        mReaction.setRegulation(regulation);
                    }
                    ModelReaction im = new ModelReaction(e, mReaction, parameters);
                    allthereactions.add(im);
                }else{
                    mReaction.setEnzyme(e);
                    ModelReaction im = new ModelReaction(e, mReaction, parameters);
                    allthereactions.add(im);
                }
            }
            
        }
        Object[][] reactioninfo = new Object[allthereactions.size()][6];
        for(int i =0;i< allthereactions.size();i++){
            reactioninfo[i][0]=i+1;
            reactioninfo[i][1]=allthereactions.get(i).getEnzyme().getName();
            if(allthereactions.get(i).getRegulation()>1){
                switch (allthereactions.get(i).getRegulation()) {
                    case 2:
                        reactioninfo[i][2]=allthereactions.get(i).getModifier().get(0).getName();
                        reactioninfo[i][3]="N/A";
                        break;
                    case 3:
                        reactioninfo[i][2]="N/A";
                        reactioninfo[i][3]=allthereactions.get(i).getModifier().get(0).getName();
                        break;
                    case 4:
                        reactioninfo[i][2]=allthereactions.get(i).getModifier().get(0).getName();
                        reactioninfo[i][3]=allthereactions.get(i).getModifier().get(1).getName();
                        break;
                    default:
                        break;
                }
            }else {
                reactioninfo[i][2]="N/A";
                reactioninfo[i][3]="N/A";
            }
            
            reactioninfo[i][4]=allthereactions.get(i).toString();
            reactioninfo[i][5]=allthereactions.get(i).getPkinetics();
            
        }
        
    }
    
    private Enzyme oldSBML(MetabolicReaction mReaction, Reaction reaction, ListOf modifiers, Enzyme e){
        
        System.out.println("there is more than 1 modifier in reaction: "+ reaction.getId());    
                
        String currentreaction = mReaction.toString();
        String message = "More than 1 modifier found in this reaction. \n"+currentreaction+"\nPlease select 1 enzyme for this reaction.";

        ArrayList<String> modifierlist = new ArrayList<String>();
        ArrayList<Species> modifierlist2 = new ArrayList<Species>();
        for (int l = 0; l <modifiers.size();l++){
            ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
            Species s = msr.getSpeciesInstance();
            e = new Enzyme(s);
            modifierlist.add(e.getName());
            modifierlist2.add(s);
        }
        
        String message2 = "More than 1 modifier found in this reaction. \n"+currentreaction+"\nIs this modifier ("+ modifierlist.get(0)+") an enzyme?";
        String modifierId;
        
        if(modifierlist.size()==1){
            Object[] options = {"Yes", "No"}; 
            int input;
            input = (int)JOptionPane.showOptionDialog(null, message2,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);
            if(input==1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
            }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }

        if (modifierlist.size()==2){

           Object[] options = {modifierlist.get(0),modifierlist.get(1), "None of the them."}; 

           int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 2 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }

        }else if (modifierlist.size()==3){
            Object[] options = {modifierlist.get(0),modifierlist.get(1),modifierlist.get(2), "None of the them."}; 
            int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 3 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }else if (modifierlist.size()==4){
            Object[] options = {modifierlist.get(0),modifierlist.get(1),modifierlist.get(2),modifierlist.get(3), "None of the them."}; 
            int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 4 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }else if (modifierlist.size()==5){
            Object[] options = {modifierlist.get(0),modifierlist.get(1),modifierlist.get(2),modifierlist.get(3),modifierlist.get(4), "None of the them."}; 
            int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 5 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }
        return e;
    }
    
    private void addSubstratesProductsFromReaction(Compound com){
        if (isinCompoudList(com) == false ) {
            allthespecies.add(com);
        }
    }
    
    public boolean isinCompoudList(Compound c) {
        boolean isin = false;

        for (int i = 0; i < allthespecies.size(); i++) {
            if (c.equals(allthespecies.get(i)) || c.getID().matches(allthespecies.get(i).getID())) {
                isin = true;
            }
        }
        return isin;
    }
    
    private void addEnzymesFromReaction(Enzyme enz){
        
        if(isinEnzymeList(enz) == false ){
            alltheenzymes.add(enz);
        }
    }
    
    public boolean isinEnzymeList(Enzyme e){
        boolean isin = false;
        
        for (int i = 0; i <alltheenzymes.size(); i++){
            if (e.equals(alltheenzymes.get(i)) || e.getID().matches(alltheenzymes.get(i).getID())){
                isin =true;
            }
        }
        return isin;
    }
    
    public ListOf<Reaction> getreactions(){
        return modelReaction;
    }
    
    public ArrayList<ModelReaction> returnReactionList(){
        return allthereactions;
    }
    
    public ArrayList<Compound> returncompounds(){
        return allthespecies;
    }
    
    public ArrayList<Enzyme> returnenzymes(){
        return alltheenzymes;
    }
    
    public String returnname(){
        return modelName;        
    }
    
    public ImportSBML returnmodel(){
        return sbmlmodel;
    }
    
    public void modparameters(double[] parameters){
        
        int track =0;
        for (int i = 0; i<reactions.size();i++){
            KL=reactions.get(i).getKineticLaw();
            LLP=KL.getListOfLocalParameters();                  
            
            for (int j = 0; j< LLP.size(); j++){
                LLP.get(j).setValue(parameters[track]);
                track++;
            }
        }
    }
    
    public void exportsbml(File file) throws XMLStreamException, SBMLException, IOException{
        SBMLWriter.write(sbmldoc, file, "output", "X");
    }
    
    public String[] getParameterNames(){
        parameternames.clear();
        for(int i=0; i < reactions.size(); i++){
            KineticLaw kl = reactions.get(i).getKineticLaw();
            ListOf lp = kl.getListOfLocalParameters();
            
            for(int j =0;j<lp.size(); j++){
                LocalParameter para = (LocalParameter) lp.get(j).clone();
                String id = para.getId();
                parameternames.add(id);
            }
            
        }
        
        String[] array = parameternames.toArray(new String[parameternames.size()]);
        
        return array;
    }
    
    public String[] getReactionID(){
        reactionIDlist.clear();
        for(int i=0; i < reactions.size(); i++){
            KineticLaw kl = reactions.get(i).getKineticLaw();
            ListOf lp = kl.getListOfLocalParameters();
            
            for(int j =0;j<lp.size(); j++){
                LocalParameter para = (LocalParameter) lp.get(j).clone();
                String id = para.getId();
                reactionIDlist.add(reactions.get(i).getId());
            }
        }
        
        String[] array = reactionIDlist.toArray(new String[reactionIDlist.size()]);
        
        return array;
    }
    
    public String[] getReactionID2(){
        String[] rID = new String[allthereactions.size()];
        
        for(int i=0;i<rID.length;i++){
            rID[i]=allthereactions.get(i).getReactionID();
        }
        
        return rID;
    }
    
    public String[] getProteinID(){
        int pkcount = 0;
        
        for(ModelReaction reaction : allthereactions){
            if(reaction.getPkinetics()==true){
                pkcount++;
            }
        }        
        String[] proteinID = new String[pkcount];
        int tick=0;
        for(ModelReaction reaction : allthereactions){
            if(reaction.getPkinetics()==true){
                proteinID[tick]=reaction.getEnzyme().getID();
                tick++;
            }
        }
        return proteinID;
    }
}
