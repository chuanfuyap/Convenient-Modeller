/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Species;

import java.util.ArrayList;
import org.sbml.jsbml.Species;

/**
 *
 * @author chuanfuyap
 */
public class Enzyme {
    
    private static final String sID = "E";
    private Species species;
    private String name, ID;
    private double concentration;
    private boolean boundaryCondition;
    private static ArrayList<String> IDArray = new ArrayList();
    private boolean proteinkinetics=true;
    
    public Enzyme(String name, double concentration) {
        
        this.ID=generateID();
        this.name = name;
        this.concentration = concentration;
        boundaryCondition = false;
    }
    
    //used when importing sbml file and creating a compound object using a species from the sbmlmodel.
    public Enzyme(Species species) {
        if (species != null) {
            this.species = species;
            this.ID = species.getId();
            addID(species.getId());
            this.name = species.getName();
            
            if(Double.isNaN(species.getInitialConcentration()) && Double.isNaN(species.getInitialAmount())){
                this.concentration = 1.0;
            }else if (Double.isNaN(species.getInitialAmount())){
                this.concentration = species.getInitialConcentration();
            }else if (Double.isNaN(species.getInitialConcentration())){
                this.concentration = species.getInitialAmount();
            }
            
            boundaryCondition = species.isBoundaryCondition();
        }
    }
    
    private static String generateID(){
        boolean isin = false;
        int j = IDArray.size()+1;
        String tempID = sID + j;
        
        for (String arrayID : IDArray) {
            if (tempID.matches(arrayID) || tempID.equals(arrayID)) {            // || here means OR
                isin = true;
                break;
            }
        }
        
        if (!isin) {
            addID(tempID);
            return tempID;
        }else{
            return ("hello");
        }
    }
    
    private static void addID(String ID) {
        IDArray.add(ID);
    }
    
    public String getName() {
        return name;
    }
    
    public double getConcentration() {
        return concentration;
    }
    
    public String getID() {
        return ID;
    }
    
    public boolean getBoundaryCondition() {
        return boundaryCondition;
    }
    
    public Species getSpecies(){
        return species;
    }
    
    public void setName(String newname) {
        name = newname;
    }
    
    public void setID(String id){
        ID = id;
    }
    
    public void setConcentration(double newconcentration) {
        concentration = newconcentration;
    }
    
    public void setBoundaryCondition(boolean boundaryCondition) {
        this.boundaryCondition = boundaryCondition;
    }
    
    public void setProteinKinetics(boolean yesno){
        proteinkinetics=yesno;
    }
    public boolean getProteinKinetics(){
        return proteinkinetics;
    }
    
}
