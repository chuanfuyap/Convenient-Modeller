/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SBML;

import java.io.File;
import java.io.IOException;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.*;

/**
 *
 * @author chuanfuyap
 */
public class ImportSBML {
    
    private SBMLDocument sbmldoc;
    private ListOf<Species> species;    
    private ListOf<Reaction> reactions;
    private ListOf<Compartment> compartments;
    private Model model;
    private String modelName;
    
    //reads a sbml file and retrieve the model (and its name) and from the model, extract the list of the reactions, speices, and compartments. 
    public ImportSBML(File fileName) throws XMLStreamException, IOException {
        sbmldoc = new SBMLDocument (SBMLReader.read(fileName));
        model = sbmldoc.getModel();
        species = (ListOf) model.getListOfSpecies();
        reactions = (ListOf) model.getListOfReactions();
        modelName = model.getName();
        compartments = (ListOf) model.getListOfCompartments();
    }
    
    public String getModelName() {
        return modelName;
    }

    public ListOf getCompartments() {
        return compartments;
    }

    public ListOf getSpecies() {
        return species;
    }
    
    public ListOf getReactions() {
        return reactions;
    }

    public Model getModel() {
        return model;
    }

    public SBMLDocument getSBMLDocument() {
        return sbmldoc;
    }

    public String toString() {
        return sbmldoc.getModel().toString();
    }
}
