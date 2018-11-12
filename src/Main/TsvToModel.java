/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import Reaction.MetabolicReaction;
import SBML.ExportSBML;
import Species.Compound;
import Species.Enzyme;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLWriter;

/**
 *
 * @author chuanfuyap
 */
public class TsvToModel {
    
    ArrayList<Compound> metabolites = new ArrayList();
    ArrayList<Enzyme> enzymes = new ArrayList();
    ArrayList<ModelReaction>  reactions = new ArrayList();
    
    public TsvToModel(){}
    
    public TsvToModel(String filelink, String outputFileLink) throws IOException, XMLStreamException{

        Path path = FileSystems.getDefault().getPath(filelink);
        List<String> datas = Files.readAllLines(path);
        
        datas.stream().map(s -> s.trim()).filter(s -> !s.isEmpty()).forEach((line) -> {
            String[] temp = line.split("\\s");
            if(temp[0].charAt(0)!='#'){
                if(temp.length==3){
                    make_and_add_metabolites(temp);
                }
                if(temp.length==2){
                    make_and_add_enzymes(temp);
                }
                if(temp.length>3){
                    make_and_add_reactions(temp);
                }
            }            
        });
        ExportSBML modelConvertor = new ExportSBML(reactions, metabolites, enzymes, "test_model");
        SBMLDocument sbml2doc = new SBMLDocument(modelConvertor.getsbmldoc());
        SBMLWriter.write(sbml2doc, new File(outputFileLink) , "output", "X");
    }
    
    public void build_model_from_tsv(String filelink) throws IOException, XMLStreamException{
        
        Path path = FileSystems.getDefault().getPath(filelink);
        List<String> datas = Files.readAllLines(path);
        
        datas.stream().map(s -> s.trim()).filter(s -> !s.isEmpty()).forEach((line) -> {
            String[] temp = line.split("\\s");
            if(temp[0].charAt(0)!='#'){
                if(temp.length==3){
                    make_and_add_metabolites(temp);
                }
                if(temp.length==2){
                    make_and_add_enzymes(temp);
                }
                if(temp.length>3){
                    make_and_add_reactions(temp);
                }
            }            
        });
        
    }
    
    private Enzyme find_enzyme(String reaction_enzyme){
        Enzyme enzyme = null;
        
        for(int i=0; i < enzymes.size();i++){
            if (reaction_enzyme.equals(enzymes.get(i).getName())){
                enzyme=enzymes.get(i);
            }
        }
        
        return enzyme;
    }
    
    private void make_and_add_metabolites(String[] input){
        String compound_name = input[0];
        double compound_concentration = Double.parseDouble(input[1]);
        Compound new_compound = new Compound(compound_name, compound_concentration);
        if(input[2].equals("TRUE")){
            new_compound.setBoundaryCondition(true);
        }
        metabolites.add(new_compound);
    }
    
    private void make_and_add_enzymes(String[] input){
        String enzyme_name = input[0];
        double enzyme_concentration = Double.parseDouble(input[1]);
        Enzyme new_enzyme = new Enzyme(enzyme_name, enzyme_concentration);
        enzymes.add(new_enzyme);
    }
    
    private void make_and_add_reactions(String[] input){
        int regulation = 1;
        String reaction_name = input[0];        
        MetabolicReaction new_mReaction = new MetabolicReaction(reaction_name);

        String activator = input[1];
        String inhibitor = input[2];
        
        ArrayList<Compound> modifiers = new ArrayList();
        if(!activator.equals("N/A")&&!inhibitor.equals("N/A")){
            regulation=4;
            modifiers.add(find_metabolite(activator));
            modifiers.add(find_metabolite(inhibitor));
        }else{
            if(activator.equals("N/A")){
                if(inhibitor.equals("N/A")){
                    regulation=1;
                }else{
                    regulation=3;
                    modifiers.add(find_metabolite(inhibitor));
                }
            }else {
                regulation=2;
                modifiers.add(find_metabolite(activator));
            } 
        }
        
        if(regulation>1){
            new_mReaction.setModifier(modifiers);
        }
        new_mReaction.setRegulation(regulation);
        
        boolean subsORprod = true;
        for (int i = 3; i<input.length;i++){
            if(input[i].equals("=")){
                subsORprod=false;
            }else{
            if(subsORprod==true){
                boolean have_stoichio = false;
                if(input[i].length() >1 && input[i].charAt(1)=='*'){
                    have_stoichio = true;
                }
                
                Compound substrate = null;
                
                if(have_stoichio==true){
                    substrate = find_metabolite(input[i].substring(2));
                    new_mReaction.addSubstrate(substrate, Character.getNumericValue(input[i].charAt(0)));
                }else{
                    substrate = find_metabolite(input[i]);
                    new_mReaction.addSubstrate(substrate);
                }
                
            }else{
                boolean have_stoichio = false;
                if(input[i].length() >1 && input[i].charAt(1)=='*'){
                    have_stoichio = true;
                }
                Compound product = null;
                if(have_stoichio==true){
                    product = find_metabolite(input[i].substring(2));
                    new_mReaction.addProduct(product, Character.getNumericValue(input[i].charAt(0)));
                }else{
                    product = find_metabolite(input[i]);
                    new_mReaction.addProduct(product);
                }
                
            }
            }
        }
        
        Enzyme reaction_enzyme = find_enzyme(reaction_name);
        ModelReaction new_reaction = new ModelReaction(reaction_enzyme, new_mReaction);

        reactions.add(new_reaction);
    }
    
    private Compound find_metabolite(String name){
        Compound metabolite = null;
        
        for(int i=0; i < metabolites.size();i++){
            if (name.equals(metabolites.get(i).getName())){
                metabolite=metabolites.get(i);
            }
        }
        
        return metabolite;
    }
    
    public ArrayList get_metabolites(){
        return metabolites;
    }
    
    public ArrayList get_enzymes(){
        return enzymes;
    }
    
    public ArrayList get_reactions(){
        return reactions;
    }
}