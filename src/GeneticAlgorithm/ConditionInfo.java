/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import java.util.HashMap;

/**
 *
 * @author chuanfuyap
 */
public class ConditionInfo {
    private HashMap metMap = new HashMap<>();
    private HashMap fluxMap = new HashMap<>();
    private HashMap Enz_Conc_to_Update = new HashMap<>();
    private HashMap BC_true_Met_to_Update = new HashMap<>();
    private HashMap solved_Met_Map = new HashMap<>();
    private HashMap solved_Flux_Map = new HashMap<>();    
    
    private boolean good_OR_bad_solution_forModel=false;
    
    public ConditionInfo(){
    }
    
    public void store_metabolites_info(HashMap metmap){
        this.metMap = metmap;
    }
    
    public void store_BCTrue_metabolites_info(HashMap metmap_tochange){
        this.BC_true_Met_to_Update = metmap_tochange;
    }
    public void store_fluxes_info(HashMap fluxmap){
        this.fluxMap=fluxmap;
    }
    
    public void store_proteins_info(HashMap promap){
        this.Enz_Conc_to_Update=promap;
    }
    
    public void store_Solved_metmap(HashMap solved_metmap){
        this.solved_Met_Map=solved_metmap;
    }
    
    public void store_Solved_fluxmap(HashMap solved_fluxmap){
        this.solved_Flux_Map=solved_fluxmap;
    }
    
    public HashMap get_metabolites_info(){
        return metMap;
    }
    
    public HashMap get_BCTrue_metabolites_info(){
        return BC_true_Met_to_Update;
    }
    public HashMap get_fluxes_info(){
        return fluxMap;
    }
    
    public HashMap get_proteins_info(){
        return Enz_Conc_to_Update;
    }
    
    public HashMap get_Solved_metmap(){
        return solved_Met_Map;
    }
    
    public HashMap get_Solved_fluxmap(){
        return solved_Flux_Map;
    }
    public void store_goodsolution(boolean yesno){
        this.good_OR_bad_solution_forModel = yesno;
    }
    public boolean goodsolution(){
        return good_OR_bad_solution_forModel;
    }
}