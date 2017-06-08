/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SBML;

import Main.ModelReaction;
import Species.*;
import java.util.ArrayList;
import java.util.Iterator;
import org.sbml.jsbml.*;

/**
 *
 * @author chuanfuyap
 */
public class ExportSBML {
    
    private Model myModel;
    private Compartment compartment;
    private SBMLDocument sbmldoc;
    private ArrayList<ModelReaction> mReactions;
    
    //sets up SBML version to lvl2v4, sets up model with Compartment at size 1 dimension.
    //default constructor
    public ExportSBML(){
        
        try{
        sbmldoc = new SBMLDocument (2,4);       //create sbml doc
        myModel = sbmldoc.createModel();        //add model to sbml doc
        compartment =  myModel.createCompartment("Cell");   //add compartment to model within sbml doc
        compartment.setSize(1.0);           //sets size of compartment
        }catch (Exception e) {
            System.out.println("Cannot Convert Model to SBML. Error 1");
        }
    }
    
    //constructor that sets name for model
    public ExportSBML(String name) {
        this();
        try{
        
        myModel.setName(name);
        
        }
        catch (Exception e) {
            System.out.println("Cannot Convert Model to SBML. Error 2 : unable to set Name for Model");
        }
    }
    
    //constructor to add species, and modelname to the model.
    public ExportSBML(ArrayList<Compound> species, String name) {
        this(name);                     //means Model2SBMLConverter(String name)
        try{
            addSpeciesToModel(species);
        }catch (Exception e) {
            System.out.println("Cannot Convert Model to SBML. Error 3 : unable to add species to the Model");
        }
    }
    
    //constructor to add species, enzyme and modelname to the model.
    public ExportSBML(ArrayList<Compound> species, ArrayList<Enzyme> enzymes, String name) {
        this(species, name);            
        try{
            addEnzymesToModel(enzymes);
        }catch (Exception e) {
            System.out.println("Cannot Convert Model to SBML. Error 4 : unable to enzyme to the Model");
        }
    }
    
    //constructor to add reactions, species, enzymes and modelname to the model. as well as a boolean to save or not.
    public ExportSBML(ArrayList<ModelReaction> mReactions, ArrayList<Compound> species, ArrayList<Enzyme> enzymes, String name) {
        this(species, enzymes, name); 
        try{            
            this.mReactions = mReactions;
            addReactionsToModel(mReactions);
        }catch (Exception e) {
            System.out.println("Cannot Convert Model to SBML. Error 5 : unable to add reactions to the Model");
        }
    }
    public ExportSBML(ArrayList<ModelReaction> mReactions, ArrayList<Compound> species, ArrayList<Enzyme> enzymes, String name, double[] parameters) {
        this(species, enzymes, name); 
        try{            
            this.mReactions = mReactions;
            addReactionsToModel(mReactions, parameters);
        }catch (Exception e) {
            System.out.println("Cannot Convert Model to SBML. Error 5 : unable to add reactions to the Model");
        }
    }
    //method to add species created to the model.
    private void addSpeciesToModel(ArrayList<Compound> speciesList) {
        Iterator iterator = speciesList.iterator();
        while (iterator.hasNext()) {
            Compound c = (Compound) iterator.next();
            Species s = myModel.createSpecies(c.getID(), c.getName(), compartment);
            s.setBoundaryCondition(c.getBoundaryCondition());
            s.setInitialConcentration(c.getConcentration());
            s.setSBOTerm(245);                                  //set SBOTerm to 245 which is macromolecule
        }
    }
    
    //method to add enzymes to the model
    private void addEnzymesToModel(ArrayList<Enzyme> enzymeList) {
        Iterator eIterator = enzymeList.iterator();
        while (eIterator.hasNext()) {
            Enzyme c = (Enzyme) eIterator.next();
            Species s = myModel.createSpecies(c.getID(), c.getName(), compartment);
            s.setBoundaryCondition(c.getBoundaryCondition());
            s.setInitialConcentration(c.getConcentration());
            s.setSBOTerm(252);                                  //set SBOTerm to 252 SBOTerm for polypeptide chain;
        }
    }
    
    private Model addReactionsToModel(ArrayList<ModelReaction> alltheReactions) {
        
        Iterator it = alltheReactions.iterator();
        while (it.hasNext()) {
            ModelReaction modelreaction = (ModelReaction) it.next();

            //add list of reactions
            Reaction reaction = new Reaction(modelreaction.getMyReaction().getReactionID());        
            reaction.setName(modelreaction.getMyReaction().getName());
            reaction.setReversible(true);
            
            //add reactants of the reaction  
            ArrayList substratesStoic = modelreaction.getSubstratesStoichio();
            ArrayList reactants = modelreaction.getSubstrates();

            try {
                for (int i = 0; i < reactants.size(); i++) {
                    Compound compound = (Compound) reactants.get(i);
                    Species s = new Species(compound.getID());
                    s.setName(compound.getName());
                    SpeciesReference sf = new SpeciesReference(s);                

                    int stoic = Integer.parseInt(substratesStoic.get(i).toString());
                    sf.setStoichiometry(stoic);
                    reaction.addReactant(sf);
                
                }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 1 "+ reaction.getId());
            }
            
            //adding enzyme as modifier to the model. 
            try{
                    Species enzyme = new Species(modelreaction.getEnzyme().getID());
                    enzyme.setName(modelreaction.getEnzyme().getName());
                    ModifierSpeciesReference enzymemod = new ModifierSpeciesReference (enzyme);
                    reaction.addModifier(enzymemod);
                
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 2 "+modelreaction.getName());
            }

            //add the products of the reaction
            ArrayList productsStoic = modelreaction.getProductsStoichio();
            ArrayList products = modelreaction.getProducts();

            try{
                
                for (int i = 0; i < products.size(); i++) {
                    Compound compound = (Compound) products.get(i);
                    Species s = new Species(compound.getID());
                    s.setName(compound.getName());
                    SpeciesReference sf = new SpeciesReference(s);                

                    int stoic = Integer.parseInt(productsStoic.get(i).toString());
                    sf.setStoichiometry(stoic);
                    reaction.addProduct(sf);
                    
                }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 3 "+ reaction.getId());
            }

            //get the kineticLaw of the reaction
            //add the transcription and translation equations
            
            int regulation = modelreaction.getRegulation();
            
            if(regulation>1){
                if(regulation<4){
                    ArrayList<Compound> modifier = modelreaction.getModifier();
                    Species mod = new Species(modifier.get(0).getID());
                    mod.setName(mod.getName());
                    ModifierSpeciesReference modsf = new ModifierSpeciesReference(mod);
                    reaction.addModifier(modsf);
                }else if(regulation==4){
                    ArrayList<Compound> modifier = modelreaction.getModifier();
                    Species mod = new Species(modifier.get(0).getID());
                    mod.setName(mod.getName());
                    ModifierSpeciesReference modsf = new ModifierSpeciesReference(mod);
                    reaction.addModifier(modsf);
                    Species mod2 = new Species(modifier.get(1).getID());
                    mod2.setName(mod2.getName());
                    ModifierSpeciesReference modsf2 = new ModifierSpeciesReference(mod2);
                    reaction.addModifier(modsf2);
                }
            }
            
            KineticLaw kineticLaw= null;
            
            try{
                kineticLaw = modelreaction.getKineticLaw();
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 4 " + modelreaction.getName());
            }
            
            try{                
                reaction.setKineticLaw(kineticLaw);
                ArrayList<Parameter> tempparameters = modelreaction.getParameters();
            
                //for loop to convert parameters to local parameter that are subsequently added to the kineticlaw.
                for (int i=0; i< tempparameters.size(); i++){
                    LocalParameter localparameter = new LocalParameter(tempparameters.get(i));                    
                    kineticLaw.addLocalParameter(localparameter);
            }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 5 "+ reaction.getId());
                System.out.println(e);
            }
            
            try{
                myModel.addReaction(reaction);
                reaction.setSBOTerm(176);
                reaction.getKineticLaw().setSBOTerm(268);                   //268 SBO enzymatic rate law 
                for (int y = 0; y < reaction.getListOfReactants().size(); y++){
                    SpeciesReference reactant = reaction.getListOfReactants().get(y);
                    reactant.setSBOTerm(15);
                }
                for (int z = 0; z < reaction.getListOfProducts().size(); z++){
                    SpeciesReference reactant = reaction.getListOfProducts().get(z);
                    reactant.setSBOTerm(11);
                }
                reaction.getListOfModifiers().get(0).setSBOTerm(460);           //460 SBO for enzymatic catalyst
                switch (regulation) {
                    case 2:
                        reaction.getListOfModifiers().get(1).setSBOTerm(534);      //533 SBO for catalytic activator;
                        break;
                    case 3:
                        reaction.getListOfModifiers().get(1).setSBOTerm(20);        //20 SBO for inhibitor
                        break;
                    case 4:
                        reaction.getListOfModifiers().get(1).setSBOTerm(534);      //533 SBO for catalytic activator;
                        reaction.getListOfModifiers().get(2).setSBOTerm(20);        //20 SBO for inhibitor
                        break;
                    default:
                        break;
                }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 6 "+modelreaction.getName());
            }
            
//            if (modelreaction.getPkinetics()==true){
//                String proteinname = modelreaction.getEnzyme().getName().replaceAll("\\W", "").replaceAll("\\s", "");
//                Reaction genereaction = new Reaction("mRNA"+proteinname);        
//                genereaction.setName(proteinname+" transcription");
//                genereaction.setReversible(false);
//                
//                Enzyme enz = modelreaction.getEnzyme();
//                Species gene = new Species("mRNA"+enz.getID());
//                gene.setName("mRNA "+enz.getName());
//                SpeciesReference genesf = new SpeciesReference(gene);
//                genesf.setStoichiometry(1);
//                genereaction.addProduct(genesf);
//                
//                KineticLaw genekineticLaw= null;
//                genekineticLaw = modelreaction.getGeneKinetics();
//                
//                genereaction.setKineticLaw(genekineticLaw);
//                ArrayList<Parameter> tempparameters3 = modelreaction.getGKparameters();
//                
//                //for loop to convert parameters to local parameter that are subsequently added to the kineticlaw.
//                for (int i=0; i< tempparameters3.size(); i++){
//                    LocalParameter localparameter = new LocalParameter(tempparameters3.get(i));
//                    genekineticLaw.addLocalParameter(localparameter);
//                }
//                
//                Species s1 = myModel.createSpecies("mRNA"+enz.getID(), "mRNA "+enz.getName(), compartment);
//                s1.setBoundaryCondition(false);
//                s1.setInitialConcentration(1);
//                s1.setSBOTerm(250);                 //250 SBO term for ribonucleic acid
//                
//                myModel.addReaction(genereaction);
//                genereaction.setSBOTerm(342);                   //342 SBOterm for molecular or genetic interaction;
//                genereaction.getKineticLaw().setSBOTerm(44);                 //44 SBOterm for mass action rate law for irreversible reactions
//                
//            }
            if (modelreaction.getPkinetics()==true){
                String proteinname = modelreaction.getEnzyme().getName().replaceAll("\\W", "").replaceAll("\\s", "");
                
                Reaction proteinreaction = new Reaction("P"+proteinname);        
                proteinreaction.setName(proteinname+" translation");
                proteinreaction.setReversible(false);
                
                Enzyme enz = modelreaction.getEnzyme();
                
//                Species gene = new Species("mRNA"+enz.getID());
////                gene.setName("mRNA "+enz.getName());
//                ModifierSpeciesReference genesf = new ModifierSpeciesReference(gene);
                Species prot = new Species(enz.getID());
                prot.setName(enz.getName());
                SpeciesReference protsf = new SpeciesReference(prot);
                protsf.setStoichiometry(1);
                
//                proteinreaction.addModifier(genesf);
                proteinreaction.addProduct(protsf);
                
                KineticLaw proteinkineticLaw= null;
                proteinkineticLaw = modelreaction.getProteinKinetics();
                
                proteinreaction.setKineticLaw(proteinkineticLaw);
                ArrayList<Parameter> tempparameters2 = modelreaction.getPKparameters();
                
                //for loop to convert parameters to local parameter that are subsequently added to the kineticlaw.
                for (int i=0; i< tempparameters2.size(); i++){
                    LocalParameter localparameter = new LocalParameter(tempparameters2.get(i));
                    proteinkineticLaw.addLocalParameter(localparameter);
                }
                myModel.addReaction(proteinreaction);
                proteinreaction.setSBOTerm(342);                                //342 SBOterm for molecular or genetic interaction;
                proteinreaction.getKineticLaw().setSBOTerm(44);                 //44 SBOterm for mass action rate law for irreversible reactions
            }

        }
        return myModel;
    }

    private Model addReactionsToModel(ArrayList<ModelReaction> alltheReactions, double[] parameters) {
        
        int track=0;
        Iterator it = alltheReactions.iterator();
        while (it.hasNext()) {
            ModelReaction modelreaction = (ModelReaction) it.next();

            //add list of reactions
            Reaction reaction = new Reaction(modelreaction.getMyReaction().getReactionID());        
            reaction.setName(modelreaction.getMyReaction().getName());
            reaction.setReversible(true);
            
            //add reactants of the reaction  
            ArrayList substratesStoic = modelreaction.getSubstratesStoichio();
            ArrayList reactants = modelreaction.getSubstrates();

            try {
                for (int i = 0; i < reactants.size(); i++) {
                    Compound compound = (Compound) reactants.get(i);
                    Species s = new Species(compound.getID());
                    s.setName(compound.getName());
                    SpeciesReference sf = new SpeciesReference(s);                

                    int stoic = Integer.parseInt(substratesStoic.get(i).toString());
                    sf.setStoichiometry(stoic);
                    reaction.addReactant(sf);
                
                }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 1 "+ reaction.getId());
            }
            
            //adding enzyme as modifier to the model. 
            try{
                    Species enzyme = new Species(modelreaction.getEnzyme().getID());
                    enzyme.setName(modelreaction.getEnzyme().getName());
                    ModifierSpeciesReference enzymemod = new ModifierSpeciesReference (enzyme);
                    reaction.addModifier(enzymemod);
                
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 2 "+modelreaction.getName());
            }

            //add the products of the reaction
            ArrayList productsStoic = modelreaction.getProductsStoichio();
            ArrayList products = modelreaction.getProducts();

            try{
                
                for (int i = 0; i < products.size(); i++) {
                    Compound compound = (Compound) products.get(i);
                    Species s = new Species(compound.getID());
                    s.setName(compound.getName());
                    SpeciesReference sf = new SpeciesReference(s);                

                    int stoic = Integer.parseInt(productsStoic.get(i).toString());
                    sf.setStoichiometry(stoic);
                    reaction.addProduct(sf);
                    
                }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 3 "+ reaction.getId());
            }

            //get the kineticLaw of the reaction
            //add the transcription and translation equations
            
            int regulation = modelreaction.getRegulation();
            
            if(regulation>1){
                if(regulation<4){
                    ArrayList<Compound> modifier = modelreaction.getModifier();
                    Species mod = new Species(modifier.get(0).getID());
                    mod.setName(mod.getName());
                    ModifierSpeciesReference modsf = new ModifierSpeciesReference(mod);
                    reaction.addModifier(modsf);
                }else if(regulation==4){
                    ArrayList<Compound> modifier = modelreaction.getModifier();
                    Species mod = new Species(modifier.get(0).getID());
                    mod.setName(mod.getName());
                    ModifierSpeciesReference modsf = new ModifierSpeciesReference(mod);
                    reaction.addModifier(modsf);
                    Species mod2 = new Species(modifier.get(1).getID());
                    mod2.setName(mod2.getName());
                    ModifierSpeciesReference modsf2 = new ModifierSpeciesReference(mod2);
                    reaction.addModifier(modsf2);
                }
            }
            
            KineticLaw kineticLaw= null;
            
            try{
                kineticLaw = modelreaction.getKineticLaw();
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 4 " + modelreaction.getName());
            }
            
            try{                
                reaction.setKineticLaw(kineticLaw);
                ArrayList<Parameter> tempparameters = modelreaction.getParameters();
                
                //for loop to convert parameters to local parameter that are subsequently added to the kineticlaw.
                for (int i=0; i< tempparameters.size(); i++){
                    LocalParameter localparameter = new LocalParameter(tempparameters.get(i));
                    localparameter.setValue(parameters[track]);
                    System.out.println(parameters[track]);
                    track++;
                    kineticLaw.addLocalParameter(localparameter);
            }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 5 "+ reaction.getId());
                System.out.println(e);
            }
            
            try{
                myModel.addReaction(reaction);
                reaction.setSBOTerm(176);
                reaction.getKineticLaw().setSBOTerm(268);                   //268 SBO enzymatic rate law 
                for (int y = 0; y < reaction.getListOfReactants().size(); y++){
                    SpeciesReference reactant = reaction.getListOfReactants().get(y);
                    reactant.setSBOTerm(15);
                }
                for (int z = 0; z < reaction.getListOfProducts().size(); z++){
                    SpeciesReference reactant = reaction.getListOfProducts().get(z);
                    reactant.setSBOTerm(11);
                }
                reaction.getListOfModifiers().get(0).setSBOTerm(460);           //460 SBO for enzymatic catalyst
                switch (regulation) {
                    case 2:
                        reaction.getListOfModifiers().get(1).setSBOTerm(534);      //533 SBO for catalytic activator;
                        break;
                    case 3:
                        reaction.getListOfModifiers().get(1).setSBOTerm(20);        //20 SBO for inhibitor
                        break;
                    case 4:
                        reaction.getListOfModifiers().get(1).setSBOTerm(534);      //533 SBO for catalytic activator;
                        reaction.getListOfModifiers().get(2).setSBOTerm(20);        //20 SBO for inhibitor
                        break;
                    default:
                        break;
                }
            }catch (Exception e) {
                System.out.println("Cannot Build Equation. Error 6 "+modelreaction.getName());
            }
            
//            if (modelreaction.getPkinetics()==true){
//                String proteinname = modelreaction.getEnzyme().getName().replaceAll("\\W", "").replaceAll("\\s", "");
//                Reaction genereaction = new Reaction("mRNA"+proteinname);        
//                genereaction.setName(proteinname+" transcription");
//                genereaction.setReversible(false);
//                
//                Enzyme enz = modelreaction.getEnzyme();
//                Species gene = new Species("mRNA"+enz.getID());
//                gene.setName("mRNA "+enz.getName());
//                SpeciesReference genesf = new SpeciesReference(gene);
//                genesf.setStoichiometry(1);
//                genereaction.addProduct(genesf);
//                
//                KineticLaw genekineticLaw= null;
//                genekineticLaw = modelreaction.getGeneKinetics();
//                
//                genereaction.setKineticLaw(genekineticLaw);
//                ArrayList<Parameter> tempparameters3 = modelreaction.getGKparameters();
//                
//                //for loop to convert parameters to local parameter that are subsequently added to the kineticlaw.
//                for (int i=0; i< tempparameters3.size(); i++){
//                    LocalParameter localparameter = new LocalParameter(tempparameters3.get(i));
//                    genekineticLaw.addLocalParameter(localparameter);
//                }
//                
//                Species s1 = myModel.createSpecies("mRNA"+enz.getID(), "mRNA "+enz.getName(), compartment);
//                s1.setBoundaryCondition(false);
//                s1.setInitialConcentration(1);
//                s1.setSBOTerm(250);                 //250 SBO term for ribonucleic acid
//                
//                myModel.addReaction(genereaction);
//                genereaction.setSBOTerm(342);                   //342 SBOterm for molecular or genetic interaction;
//                genereaction.getKineticLaw().setSBOTerm(44);                 //44 SBOterm for mass action rate law for irreversible reactions
//                
//            }
            if (modelreaction.getPkinetics()==true){
                String proteinname = modelreaction.getEnzyme().getName().replaceAll("\\W", "").replaceAll("\\s", "");
                
                Reaction proteinreaction = new Reaction("P"+proteinname);        
                proteinreaction.setName(proteinname+" translation");
                proteinreaction.setReversible(false);
                
                Enzyme enz = modelreaction.getEnzyme();
                
//                Species gene = new Species("mRNA"+enz.getID());
////                gene.setName("mRNA "+enz.getName());
//                ModifierSpeciesReference genesf = new ModifierSpeciesReference(gene);
                Species prot = new Species(enz.getID());
                prot.setName(enz.getName());
                SpeciesReference protsf = new SpeciesReference(prot);
                protsf.setStoichiometry(1);
                
//                proteinreaction.addModifier(genesf);
                proteinreaction.addProduct(protsf);
                
                KineticLaw proteinkineticLaw= null;
                proteinkineticLaw = modelreaction.getProteinKinetics();
                
                proteinreaction.setKineticLaw(proteinkineticLaw);
                ArrayList<Parameter> tempparameters2 = modelreaction.getPKparameters();
                
                //for loop to convert parameters to local parameter that are subsequently added to the kineticlaw.
                for (int i=0; i< tempparameters2.size(); i++){
                    LocalParameter localparameter = new LocalParameter(tempparameters2.get(i));
                    proteinkineticLaw.addLocalParameter(localparameter);
                }
                myModel.addReaction(proteinreaction);
                proteinreaction.setSBOTerm(342);                                //342 SBOterm for molecular or genetic interaction;
                proteinreaction.getKineticLaw().setSBOTerm(44);                 //44 SBOterm for mass action rate law for irreversible reactions
            }

        }
        return myModel;
    }
    public Model getModel() {
        return myModel;
    }
    
    public SBMLDocument getsbmldoc(){
        return sbmldoc;
        
    }
    
    public void modparameters(double[] parameters){
        
        int track =0;
        ListOf<LocalParameter> LLP;
        for (int i = 0; i<myModel.getListOfReactions().size();i++){
            KineticLaw KL=myModel.getListOfReactions().get(i).getKineticLaw();
            LLP=KL.getListOfLocalParameters();                  
            
            for (LocalParameter LLP1 : LLP) {
                LLP1.setValue(parameters[track]);
                track++;
            }
        }
    }
    
}
