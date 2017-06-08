/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SteadyState;

import GeneticAlgorithm.ReadData;
import GeneticAlgorithm.SystemToSolve;
import Main.ModelReaction;
import Main.SBMLFile;
import Species.Compound;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.ListOf;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class MultiSS {
    
    double[][] parameter_matrix;
    double[][] metabolite_matrix;
    double[][] fluxes_matrix;
    int parameter_row;
    int parameter_col;
    int metabolite_row;
    int metabolite_col;
    int fluxes_row;
    int fluxes_col;
    double[] fitness_array;
    double[] rmse_array;
    
    public MultiSS (String modellink, String inputlink,  int start, int end) throws XMLStreamException, IOException, ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException{
        int range= (end - start);
        parameter_col=range;
        metabolite_col=range;
        fluxes_col=range;
        fitness_array = new double[range];
        rmse_array = new double[range];
        
        String link = modellink+start+".xml";
        ReadData data = new ReadData(new File(inputlink));
        SystemToSolve sys = new SystemToSolve(link, data);
        File SBMLFile = new File(link);
        SBMLFile modelreactions = new SBMLFile(SBMLFile);

        
        ArrayList<Compound> temp_Compounds = modelreactions.returncompounds();
        int tick=0;
        for(Compound comp : temp_Compounds){
            if (comp.getBoundaryCondition()==false){
                tick++;
            }
        }
        
        parameter_row = sys.getParametersCount();
        metabolite_row = tick;
        fluxes_row = modelreactions.returnReactionList().size();
        
        parameter_matrix = new double[parameter_row][parameter_col];
        metabolite_matrix = new double[metabolite_row][metabolite_col];
        fluxes_matrix = new double[fluxes_row][fluxes_col];
        
        int count=0;
        for(int i=start; i<end; i++){
            link = modellink+i+".xml";
            SBMLFile = new File(link);
            modelreactions = new SBMLFile(SBMLFile);
            ArrayList<Compound> Compounds = modelreactions.returncompounds();
            ArrayList<ModelReaction> Reactions = modelreactions.returnReactionList();
            
            
            sys = new SystemToSolve(link, data);
            
            double[] parameter = new double[sys.getParametersCount()];
            
            
            ListOf<Reaction> reactionlist = modelreactions.getreactions();
            int counter=0;
            KineticLaw kl;
            for(int j = 0; j<reactionlist.size();j++){

                kl = reactionlist.get(j).getKineticLaw();

                for(int k = 0; k<kl.getListOfLocalParameters().size();k++){
                   parameter[counter]=kl.getListOfLocalParameters().get(k).getValue();
                   String name = kl.getListOfLocalParameters().get(k).getId();
                   counter++;
                }
            }
            
            for(int k = 0; k<sys.getParametersCount();k++){
                parameter_matrix[k][count]=parameter[k];
            }
            
            tick=0;
            for(Compound comp : Compounds){
                if (comp.getBoundaryCondition()==false){
                    tick++;
                }
            }

            double[] activeconc = new double[tick];
            int tock=0;
            for(Compound comp : Compounds){
                if (comp.getBoundaryCondition()==false){
                    activeconc[tock]=comp.getConcentration();
                    tock++;
                }
            }
            
            sys.runEst(parameter);
            double[] metabolites=sys.get_MetVector();
            
            for(int k = 0; k<metabolites.length;k++){
                metabolite_matrix[k][count]=metabolites[k];
            }
            
            double[] fluxes= new double[fluxes_row];
            
            HashMap map = sys.getEstFluxMap();
            int fluxcount=0;
            for(ModelReaction react : Reactions ){
                fluxes[fluxcount]=(double) map.get(react.getMyReaction().getReactionID());
                fluxcount++;
            }
            
            for(int k = 0; k<fluxes.length;k++){
                fluxes_matrix[k][count]=fluxes[k];
            }
            
            fitness_array[count]=get_fitness(sys);
            rmse_array[count]=get_rmse(sys);
            count++;
        }
    }
    
    public double[][] get_parameter_matrix(){
        return parameter_matrix;
    }
    
    public int get_para_row(){
        return parameter_row;
    }
    
    public int get_para_col(){
        return parameter_col;
    }
    
    public double[][] get_metabolite_matrix(){
        return metabolite_matrix;
    }
    
    public int get_metabolite_row(){
        return metabolite_row;
    }
    
    public int get_metabolite_col(){
        return metabolite_col;
    }
    
    public double[][] get_fluxes_matrix(){
        return fluxes_matrix;
    }
    
    public int get_fluxes_row(){
        return fluxes_row;
    }
    
    public int get_fluxes_col(){
        return fluxes_col;
    }
    
    private double get_fitness(SystemToSolve system){
        double summed_error = 0.0;
        double f;
        
        HashMap estMet = system.getEstMetMap();      //estimated metconc values
        HashMap estFlux = system.getEstFluxMap();    //estimated flux values
        HashMap expMet = system.getMetMap();      //estimated metconc values
        HashMap expFlux = system.getFluxMap();    //estimated flux values
        
        ArrayList<Compound> Compounds=system.getCompounds();
        ArrayList<ModelReaction> Reactions=system.getReactions();
        
        int mettrack=0;
        int fluxtrack=0;
        
        double mettotal=0;
        double fluxtotal=0;
        
        if(expMet!=null){
            for (Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    String ID = comp.getID();
                    if(expMet.containsKey(ID)==true && (double) expMet.get(ID)>=0){
                        double difference =(double) expMet.get(ID)-(double)estMet.get(ID);
                        double percentagedifference=0;
                        if((double) expMet.get(ID)>0){
                            percentagedifference = (Math.abs(difference)/(double) expMet.get(ID))*100;
                        }else{
                            percentagedifference = Math.abs(difference)*100;
                        }                        
                        mettotal+=Math.pow(percentagedifference, 2);
                        mettrack++;
                    }
                }
            }
            }
            
            if(expFlux!=null){
            for(ModelReaction reaction : Reactions){
                String ID = reaction.getReactionID();
                if(expFlux.containsKey(ID)==true && (double) expFlux.get(ID)>=0){
                    double difference = (double) expFlux.get(ID) - (double) estFlux.get(ID);
                    double percentagedifference=0;
                    if((double) expFlux.get(ID)>0){
                        percentagedifference = (Math.abs(difference)/(double) expFlux.get(ID))*100;
                    }else{
                        percentagedifference = Math.abs(difference)*100;
                    }
                    fluxtotal+=Math.pow(percentagedifference, 2);
                    fluxtrack++;
                }
            }
            }
            if(mettrack==0){
                double fluxdifference = fluxtotal/fluxtrack;
                summed_error=fluxdifference;
            }else if (fluxtrack==0){
                double metdifference = mettotal/mettrack;
                summed_error=metdifference;
            }else{
                double metdifference = mettotal/mettrack;
                metdifference*=0.5;
                double fluxdifference = fluxtotal/fluxtrack;
                fluxdifference*=0.5;

                summed_error+=metdifference+fluxdifference;
            }

        
        f = 1.0 / (1.0 + summed_error);
        
        return f;
    }
    
    private double get_rmse(SystemToSolve system){
        double summed_error = 0.0;
        
        HashMap estMet = system.getEstMetMap();      //estimated metconc values
        HashMap estFlux = system.getEstFluxMap();    //estimated flux values
        HashMap expMet = system.getMetMap();      //estimated metconc values
        HashMap expFlux = system.getFluxMap();    //estimated flux values
        
        ArrayList<Compound> Compounds=system.getCompounds();
        ArrayList<ModelReaction> Reactions=system.getReactions();
        
        int mettrack=0;
        int fluxtrack=0;
        
        double mettotal=0;
        double fluxtotal=0;
        
        if(expMet!=null){
            for (Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    String ID = comp.getID();
                    if(expMet.containsKey(ID)==true && (double) expMet.get(ID)>=0){
                        double difference = (double) expMet.get(ID)-(double)estMet.get(ID);
                                       
                        mettotal+=Math.pow(difference, 2);
                        mettrack++;
                    }
                }
            }
            }
            
            if(expFlux!=null){
            for(ModelReaction reaction : Reactions){
                String ID = reaction.getReactionID();
                if(expFlux.containsKey(ID)==true && (double) expFlux.get(ID)>=0){
                    double difference = (double) expFlux.get(ID) - (double) estFlux.get(ID);
                   
                    fluxtotal+=Math.pow(difference, 2);
                    fluxtrack++;
                }
            }
            }
            
            if(mettrack==0){
                double fluxdifference = fluxtotal/fluxtrack;
                summed_error=fluxdifference;
            }else if (fluxtrack==0){
                double metdifference = mettotal/mettrack;
                summed_error=metdifference;
            }else{
                
                double totaldiff = mettotal + fluxtotal;
                double total = mettrack + fluxtrack;

                summed_error+=totaldiff/total;
            }
    
        return Math.sqrt(summed_error);
    }
    
    public double[] get_fit_array(){
        return fitness_array;
    }
    
    public double[] get_rmse_array(){
        return rmse_array;
    }
}
