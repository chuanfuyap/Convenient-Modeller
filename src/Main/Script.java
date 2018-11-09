/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import GeneticAlgorithm.GAlgorithm;
import GeneticAlgorithm.ReadData;
import GeneticAlgorithm.SystemToSolve;
import Species.Compound;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import odeBridge.TCodesolve;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.ListOf;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class Script {
    
    private ArrayList<Compound> Compounds;
    private ArrayList<ModelReaction> Reactions;   
    private SystemToSolve sys;
    private SBMLFile modelreactions;
    private String ModelLink;
    private ReadData data;
    
    public Script() {
        
    }
    
    public void addModel(String ModelLink) throws IOException, XMLStreamException, ModelOverdeterminedException{
        File SBMLFile = new File(ModelLink);
        this.ModelLink=ModelLink;
        
        modelreactions = new SBMLFile(SBMLFile);
        Compounds = modelreactions.returncompounds();
        Reactions = modelreactions.returnReactionList();
        
    }
    
    public void addData(String DataLink) throws IOException{
        File InputData = new File(DataLink);
        data = new ReadData(InputData);
    }
    
    public void build()  throws IOException, XMLStreamException, ModelOverdeterminedException{
        sys = new SystemToSolve(ModelLink, data);
    }
    
    public void runGA(String outputLink, int popsize, int maxgen, int plateau, int plague, double[][] Grange)throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException, XMLStreamException, IOException{
        GAlgorithm ga = new GAlgorithm(sys, true, popsize, maxgen, plateau, plague, Grange);
        ga.run();
        ArrayList params = ga.getParameters();
        ArrayList fitnessTrack = ga.getFitnessTrack();

        double[] parameters=(double[]) params.get(0);
        modelreactions.modparameters(parameters);
        File outputfile = new File (outputLink+".xml");
        modelreactions.exportsbml(outputfile);
        PrintWriter writer = new PrintWriter( outputLink+"Fitness.txt", "UTF-8" );

        for(Object fitness:fitnessTrack){
            writer.println(fitness);
        }

        writer.close();
    }
    
    public void SSsolution() throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, XMLStreamException, IOException {
        double[] parameter = new double[sys.getParametersCount()];
        ListOf<Reaction> reactionlist = modelreactions.getreactions();
        int counter=0;
        KineticLaw kl=null;
        for(int j = 0; j<reactionlist.size();j++){
            
            kl = reactionlist.get(j).getKineticLaw();
            
            for(int k = 0; k<kl.getListOfLocalParameters().size();k++){
               parameter[counter]=kl.getListOfLocalParameters().get(k).getValue();
               String name = kl.getListOfLocalParameters().get(k).getId();
               counter++;
            }
        }
        System.out.println("Parameter Count:\t"+counter);
        sys.solve_ODE_model(parameter,Integer.toHexString(this.hashCode()));
        if(sys.goodsolution()){
            HashMap mmap = sys.get_Condition_List().get(0).get_Solved_metmap();
            for (Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                String ID = comp.getID();

                System.out.println(comp.getName()+"\t"+ID+"\t"+mmap.get(ID));
                }
            }
            System.out.println();
            HashMap map = sys.get_Condition_List().get(0).get_Solved_fluxmap();
            for(ModelReaction react : Reactions ){
                String ID = react.getMyReaction().getReactionID();
                System.out.println(react.getMyReaction().getName()+"\t"+ID+"\t"+map.get(ID));
            }

            System.out.println("\nFitness\t"+get_fitness(sys));
            System.out.println("\nRMSE\t"+get_rmse(sys));
        }else{
            System.out.println("parameters are no good");
        }
    }
    private double get_rmse(SystemToSolve system){
        double summed_error = 0.0;
        
        HashMap estMet = system.get_Condition_List().get(0).get_Solved_metmap();      //estimated metconc values
        HashMap estFlux = system.get_Condition_List().get(0).get_Solved_fluxmap();    //estimated flux values
        HashMap expMet = system.get_Condition_List().get(0).get_metabolites_info();      //estimated metconc values
        HashMap expFlux = system.get_Condition_List().get(0).get_fluxes_info();    //estimated flux values
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
    
    private double get_fitness(SystemToSolve system){
        double summed_error = 0.0;
        double f;
        
        HashMap estMet = system.get_Condition_List().get(0).get_Solved_metmap();      //estimated metconc values
        HashMap estFlux = system.get_Condition_List().get(0).get_Solved_fluxmap();    //estimated flux values
        HashMap expMet = system.get_Condition_List().get(0).get_metabolites_info();      //estimated metconc values
        HashMap expFlux = system.get_Condition_List().get(0).get_fluxes_info();    //estimated flux values
        
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
    
    public void runTimeSeries()throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException, XMLStreamException, IOException{
        int activeCompNum=0;
        for (Compound comp: Compounds){
            if (comp.getBoundaryCondition()==false){
                activeCompNum++;
            }
        }
        String[] variables = new String[activeCompNum];
        int count=0;
            for(Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    variables[count]=comp.getID();
                    count++;
                }
            }
        TCodesolve ode = new TCodesolve(ModelLink, variables, sys.getParametersCount(), modelreactions.getParameterNames(), modelreactions.getReactionID(), modelreactions.getReactionID2());
        double[] parameter = new double[sys.getParametersCount()];
        ListOf<Reaction> reactionlist = modelreactions.getreactions();
        int counter=0;
        KineticLaw kl=null;
        for(int j = 0; j<reactionlist.size();j++){
            
            kl = reactionlist.get(j).getKineticLaw();
            
            for(int k = 0; k<kl.getListOfLocalParameters().size();k++){
               parameter[counter]=kl.getListOfLocalParameters().get(k).getValue();
               String name = kl.getListOfLocalParameters().get(k).getId();
               counter++;
            }
        }

        sys.solve_ODE_model(parameter,Integer.toHexString(this.hashCode()));
        HashMap metmapTC = sys.getTCestmetMap();
        double[] metTC = (double[])metmapTC.get("S8");
        for(double num : metTC){
            System.out.println(num);
        }
    }
    
  
    public void txtToModel(String fileLink, String outputLink, int start){
        File textfile = new File(fileLink);
        try {
            BufferedReader input = new BufferedReader (
                                        new InputStreamReader(
                                            new FileInputStream(textfile)));
            String line;
            while ((line = input.readLine()) != null) {
                String[] token = line.split("\\s");
                if(token[0].equals("KP")){
                    double[] parameters = new double[token.length-1];
                    for(int i=0; i<parameters.length;i++){
                        parameters[i]=Double.parseDouble(token[i+1]);
                    }
                    modelreactions.modparameters(parameters);
                    File outputfile = new File (outputLink+start+".xml");
                    modelreactions.exportsbml(outputfile);
                    start++;
                }
            }
            
        }catch (Exception e) {
        }
    }
    
}