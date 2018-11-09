/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author chuanfuyap
 */
public class ReadData {
    private String[] IDs;
    private HashMap TCdata = new HashMap<>(); //time course data
    private boolean SSorTC = true;      //true for steady state, false for Time course data
    private ArrayList time=new ArrayList<>();
    private ArrayList<HashMap> conditions = new ArrayList<>();
    
    //Boolean SS to determine data being estimated is of SteadyState or of Time Course data, TRUE IF SS
    public ReadData(File datafile) throws FileNotFoundException, IOException{
        
        try {
            
            BufferedReader input = new BufferedReader (
                                        new InputStreamReader(
                                            new FileInputStream(datafile)));
            
            ArrayList lines = new ArrayList();            
            
            input.lines().forEach((line) -> {
                if(line.trim().length()>0){
                    lines.add(line);
                }
            });
            
            String firstline = (String) lines.get(0);
            String[] temp = firstline.split("\\s");
            if(temp[0].equals("Time")){
                SSorTC=false;   //false for time course
            }
            
            if (SSorTC==true){  //true for steady state data
                IDs= firstline.split("\\s");
                    int conditions_count = IDs.length - 1;
                    for(int i = 0; i<conditions_count; i++){
                        HashMap ss_condition = new HashMap<>();
                        for(Object row : lines){
                            String holder = (String) row;
                            String[] token = holder.split("\\s");
                            if(!token[i+1].equals("N/A")){
                                ss_condition.put(token[0], Double.parseDouble(token[i+1]));
                            }
                        }
                        conditions.add(ss_condition);
                    }
                                    
            }else{      //handling time course data.
                IDs=firstline.split("\\s");
                int row = lines.size();
                int col = IDs.length;
                                
                String[][] timeCoursedata = new String[row][col];
                
                for (int i =0 ; i<row; i++){
                    String newline = (String) lines.get(i);
                    timeCoursedata[i] = newline.split("\\s");
                    if(i>0){
                        time.add(timeCoursedata[i][0]);
                    }
                }
                for (int i =1 ; i<col; i++){
                    double[] timeCourse = new double[row-1];
                    
                    for(int j=1; j<row;j++){
                        if(timeCoursedata[j][i].equals("N/A")||timeCoursedata[j][i].equals("null")|| timeCoursedata[j][i]==null){
                            timeCourse[j-1]=-1;
                        }else{
                            timeCourse[j-1] = Double.parseDouble(timeCoursedata[j][i]);
                        }
                    }                 
                    TCdata.put(IDs[i], timeCourse);
                }                
            }
            
        } catch (Exception e) {
            System.out.println(e);
            System.exit(1);
        }
        
    }
    public String[] getID(){
        return IDs;
    }
    
    public HashMap getTCdata(){
        return TCdata;
    }
    
    public boolean SSorTC(){
        return SSorTC;
    }
    
    public double[] returnTime(){
        double[] timecourse = new double[time.size()];
        for(int i =0; i <time.size();i++){
            timecourse[i]=Double.parseDouble((String) time.get(i));
        }
        
        return timecourse;
    }
    
    public ArrayList getMCdata(){
        return conditions;
    }
}

