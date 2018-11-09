/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import Main.ModelReaction;
import Main.SBMLFile;
import Species.Compound;
import Species.Enzyme;
import SteadyState.sosLIBssSolver;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import odeBridge.TCodesolve;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class SystemToSolve {
    private int totalParameters=0;
    private ArrayList<Compound> Compounds;
    private ArrayList<ModelReaction> Reactions;
    private ArrayList<Enzyme> Enzymes;
    private ArrayList<Integer> fRateCIndex = new ArrayList<>(); //rate constant index
    private ArrayList<Integer> rRateCIndex = new ArrayList<>();
    private ArrayList<Integer> kaRateCIndex = new ArrayList<>(); //rate constant index
    private ArrayList<Integer> kiRateCIndex = new ArrayList<>();
    private HashMap TCmetMap = new HashMap<>();
    private HashMap TCprotMap = new HashMap<>();
    private HashMap TCfluxMap = new HashMap<>();
    private HashMap TCestMetMap = new HashMap<>();
    private HashMap TCestProtMap = new HashMap<>();
    private HashMap TCestFluxMap = new HashMap<>();
    private SBMLFile modelreactions;
    private String link;
    private String originallink;
    private boolean SteadyState_OR_TimeCourse_data;
    private ReadData data;
    private boolean good_OR_bad_solution_forModel=false;
    private ArrayList<ConditionInfo> conditions_list = new ArrayList();    
    private ArrayList multi_condition_data;    
    private double[] metabolites;
    
    public SystemToSolve(){}
    
    public SystemToSolve(String Link, ReadData data) throws IOException, XMLStreamException, ModelOverdeterminedException{
        this.data=data;
        SteadyState_OR_TimeCourse_data = data.SSorTC();
        File SBMLFile = new File(Link);
        this.modelreactions = new SBMLFile(SBMLFile);
        this.Compounds=modelreactions.returncompounds();
        this.Reactions=modelreactions.returnReactionList();
        this.Enzymes=modelreactions.returnenzymes();
        this.link=Link;
        this.originallink=Link;
        
        store_parameter_count_and_make_parameterIndex();
        
        HashMap tcdata;
        if(SteadyState_OR_TimeCourse_data==true){ // for steady state data, which is most data

            multi_condition_data = data.getMCdata();

            for (int i =0;i<multi_condition_data.size();i++){
                HashMap condition_data = (HashMap) multi_condition_data.get(i);
                ConditionInfo cond_info = new ConditionInfo();
                HashMap condition_met_info = new HashMap();
                HashMap condition_bctruemet_info = new HashMap();
                HashMap condition_prot_info = new HashMap();
                HashMap condition_flux_info = new HashMap();

                boolean do_we_have_metabolites = is_there_metabolites(condition_data);
                boolean do_we_have_fluxes = is_there_fluxes(condition_data);
                boolean do_we_have_proteins = is_there_proteins(condition_data);

                //storing the initial concentration of the active compounds and the SS value from input data for PE;
                if (do_we_have_metabolites!=false){
                    for (Compound compound : Compounds){
                        if (compound.getBoundaryCondition()==false){
                            if(condition_data.containsKey(compound.getID())==true){
                                condition_met_info.put(compound.getID(), condition_data.get(compound.getID()));
                            }
                        }
                        if (compound.getBoundaryCondition()==true){
                            if(condition_data.containsKey(compound.getID())==true){
                                condition_bctruemet_info.put(compound.getID(), condition_data.get(compound.getID()));
                            }
                        }
                    }
                }
                if (do_we_have_proteins!=false){
                    for (Enzyme enz : Enzymes){
                        if(condition_data.containsKey(enz.getID())==true){
                            condition_prot_info.put(enz.getID(),condition_data.get(enz.getID()));
                        }
                    }
                }

                if (do_we_have_fluxes!=false){
                    for (ModelReaction reaction : Reactions){
                        if(condition_data.containsKey(reaction.getReactionID())==true){
                            condition_flux_info.put(reaction.getReactionID(),condition_data.get(reaction.getReactionID()));
                        }
                    }
                }else{
                    condition_flux_info =null;
                }

                cond_info.store_BCTrue_metabolites_info(condition_bctruemet_info);
                cond_info.store_metabolites_info(condition_met_info);
                cond_info.store_proteins_info(condition_prot_info);
                cond_info.store_fluxes_info(condition_flux_info);

                conditions_list.add(cond_info);
            }
        }
        
        else{   //for time course data which is severely lacking
            tcdata= data.getTCdata();
            boolean do_we_have_metabolites = is_there_metabolites(tcdata);
            boolean do_we_have_fluxes = is_there_fluxes(tcdata);
            boolean do_we_have_proteins = is_there_proteins(tcdata);
            if (do_we_have_metabolites!=false){
            for (Compound compound : Compounds){
                if(compound.getBoundaryCondition()==false){
                if(tcdata.containsKey(compound.getID())==true){
                TCmetMap.put(compound.getID(), (double[]) tcdata.get(compound.getID()));
                }
                }
            }
            }else{
                TCmetMap =null;
            }
            
            if (do_we_have_proteins!=false){
            for (Enzyme enzyme : Enzymes){
                if(tcdata.containsKey(enzyme.getID())==true){
                TCprotMap.put(enzyme.getID(), (double[]) tcdata.get(enzyme.getID()));
                }
            }
            }else{
                TCprotMap =null;
            }
                        
            if (do_we_have_fluxes!=false){
            for (ModelReaction reaction : Reactions){
                if(tcdata.containsKey(reaction.getReactionID())==true){
                    TCfluxMap.put(reaction.getReactionID(), (double[]) tcdata.get(reaction.getReactionID()));
                }
            }
            }else{
                TCfluxMap =null;
            }
        }
    }
    
    public void solve_ODE_model(double[] parameters, String mem_address) throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, XMLStreamException, IOException {
    
        if(SteadyState_OR_TimeCourse_data==true){       //solve model for steady state values
            for (int i =0; i<conditions_list.size(); i++){
                HashMap Enz_Conc_to_Update = conditions_list.get(i).get_proteins_info();
                HashMap BC_true_Met_to_Update = conditions_list.get(i).get_BCTrue_metabolites_info();
                
                if(BC_true_Met_to_Update.size()>0){
                    modelreactions.modExtMet(BC_true_Met_to_Update);
                    for (Compound compound : Compounds){
                        if (compound.getBoundaryCondition()==true){
                            if(BC_true_Met_to_Update.containsKey(compound.getID())==true){
                                compound.setConcentration((double) BC_true_Met_to_Update.get(compound.getID()));
                            }
                        }
                    }
                }
                if(Enz_Conc_to_Update.size()>0){
                    modelreactions.modEnzymes(Enz_Conc_to_Update);
                    for(ModelReaction rxn : Reactions){
                        if(Enz_Conc_to_Update.containsKey(rxn.getEnzyme().getID())==true){
                            rxn.getEnzyme().setConcentration((double) Enz_Conc_to_Update.get(rxn.getEnzyme().getID()));
                        }
                    }
                }
                
                String tmpname = mem_address;
                tmpname+="_"+i+".xml";

                File tmpxml = new File (tmpname);
                if(tmpxml.exists()){
                    String rando = Double.toHexString(RN(0.01, 10) * RN(0.01,10)).replaceAll("\\.", "_");
                    tmpname+=rando;
                    tmpxml = new File(tmpname);
                }
                
                if(BC_true_Met_to_Update.size()>0|| Enz_Conc_to_Update.size()>0){
                    modelreactions.exportsbml(tmpxml);
                    link = tmpname;
                }else{
                    link = originallink;
                }
                
                sosLIBssSolver steadystateresults = new sosLIBssSolver(Compounds, Reactions, parameters, link);

                if(steadystateresults.solveSS()==true){
                    conditions_list.get(i).store_Solved_metmap(steadystateresults.getSolvedMMap());
                    conditions_list.get(i).store_Solved_fluxmap(steadystateresults.getSolvedFMap());
                    good_OR_bad_solution_forModel=true;
                }else{
                    good_OR_bad_solution_forModel=false;
                }
                tmpxml.delete();

                    
            }
            
        }
        else{       //solve model for time series output/values
            double[] time = data.returnTime();
            double[] parameter = parameters;  
            int activeCompNum=0;
            for (Compound comp: Compounds){
                if (comp.getBoundaryCondition()==false){
                    activeCompNum++;
                }
            }

            String[] paranames = modelreactions.getParameterNames();
            String[] reactionID = modelreactions.getReactionID();
            String[] variables= new String[activeCompNum];
            String[] reactionID2 = modelreactions.getReactionID2();
            String[] proteinID = modelreactions.getProteinID();

            int count=0;
            for(Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    variables[count]=comp.getID();
                    count++;
                }
            }

            TCodesolve ode = new TCodesolve(link, variables, getParametersCount(), paranames, reactionID, reactionID2);

            double output[][] = ode.runsolver(parameter, time);

            for(int j=1;j<variables.length+1+reactionID2.length;j++){
                double[] temparray = new double[time.length];
                for(int k =0;k<time.length;k++){                   
                    temparray[k]=output[k][j];
                }
                if((j-1)<activeCompNum){
                    TCestMetMap.put(variables[j-1], temparray);
                }else{
                    TCestFluxMap.put(reactionID2[j-activeCompNum-1],temparray);
                }
            }
        }
    }
    
    public HashMap getTCmetMap(){
        return TCmetMap;
    }
    
    public HashMap getTCprotMap(){
        return TCprotMap;
    }
    
    public HashMap getTCfluxMap(){
        return TCfluxMap;
    }
    
    public HashMap getTCestmetMap(){
        return TCestMetMap;
    }
    
    public HashMap getTCestprotMap(){
        return TCestProtMap;
    }
    
    public HashMap getTCestfluxMap(){
        return TCestFluxMap;
    }
    public int[] givefRIndex(){
        int[] rateindex = new int[fRateCIndex.size()];
        for(int i=0;i<rateindex.length; i++){
            rateindex[i]=fRateCIndex.get(i);
        }
        
        return rateindex;
    }
    
    public int[] giverRIndex(){
        int[] rateindex = new int[rRateCIndex.size()];
        for(int i=0;i<rateindex.length; i++){
            rateindex[i]=rRateCIndex.get(i);
        }
        
        return rateindex;
    }
    
    public int[] giverKiIndex(){
        int[] rateindex = new int[kiRateCIndex.size()];
        for(int i=0;i<rateindex.length; i++){
            rateindex[i]=kiRateCIndex.get(i);
        }
        
        return rateindex;
    }
    
    public int[] giverKaIndex(){
        int[] rateindex = new int[kaRateCIndex.size()];
        for(int i=0;i<rateindex.length; i++){
            rateindex[i]=kaRateCIndex.get(i);
        }
        
        return rateindex;
    }
    
    public int getParametersCount(){
        return totalParameters;
    }
     
    public ArrayList<Compound> getCompounds(){
        return Compounds;
    }
    
    public ArrayList<ModelReaction> getReactions(){
        return Reactions;
    }
    
    public boolean goodsolution(){
        return good_OR_bad_solution_forModel;
    }
    
    public int timelength(){
        return data.returnTime().length;
    }
    
    private boolean is_there_metabolites(HashMap data){
        boolean is_it_there = false;
        int count=0;
        for (Compound compound : Compounds){
            if(compound.getBoundaryCondition()==false){
                String compound_ID = compound.getID();
                if(data.containsKey(compound_ID)){
                    count++;
                }
            }
        }
        if(count>0){
            is_it_there=true;
        }
        return is_it_there;
    }
    
    private boolean is_there_fluxes(HashMap data){
        boolean is_it_there = false;
        int count=0;
        for (ModelReaction reaction : Reactions){
            String reaction_ID = reaction.getReactionID();
            if(data.containsKey(reaction_ID)){
                count++;
            }
        }
        if(count>0){
            is_it_there=true;
        }
        return is_it_there;
    }
    
    private boolean is_there_proteins(HashMap data) {
        boolean is_it_there = false;
        int count=0;
        for (Enzyme enzyme : Enzymes){
            String enzyme_ID = enzyme.getID();
            if(data.containsKey(enzyme_ID)){
                count++;
            }
        }
        if(count>0){
            is_it_there=true;
        }
        return is_it_there; 
    }
    
    public boolean isitSSorTC(){
        return SteadyState_OR_TimeCourse_data;
    }
    
    public double[] get_MetVector(){
        return metabolites;
    }
    
    private void store_parameter_count_and_make_parameterIndex(){
        for(ModelReaction reaction : Reactions){
            int numSub = reaction.getSubstrates().size();
            int numProd = reaction.getProducts().size();

            int regulation = reaction.getRegulation();
            int para =2;
            if(regulation>1){
                if(regulation<4){
                    para+=1;
                }else if(regulation==4){
                    para+=2;
                }
            }
            para+=numSub;
            para+=numProd;
            fRateCIndex.add(totalParameters);
            rRateCIndex.add(totalParameters+1);
            if(regulation>1){
                switch (regulation) {
                    case 2:                            
                        kaRateCIndex.add(totalParameters+2);
                        break;
                    case 3:
                        kiRateCIndex.add(totalParameters+2);
                        break;
                    case 4:
                        kaRateCIndex.add(totalParameters+2);
                        kiRateCIndex.add(totalParameters+3);
                        break;
                    default:
                        break;
                }
            }
            totalParameters+=para;
        }
    }
    
    public ArrayList<ConditionInfo> get_Condition_List(){
        return conditions_list;
    }
    
    private double RN(double lower, double upper) {
        RandomDataGenerator random = new RandomDataGenerator();
        return random.nextUniform(lower, upper);
    }
}