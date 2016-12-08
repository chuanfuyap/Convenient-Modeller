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
    private double[] SSvalues;
    private ArrayList<double[]> TimeCourseValues= new ArrayList<>();
    private HashMap SSdata = new HashMap<>(); //steady state data
    private HashMap TCdata = new HashMap<String, double[]>(); //time course data
    private boolean SSorTC = true;
    private ArrayList time=new ArrayList<>();
    
    //Boolean SS to determine data being estimated is of SteadyState or of Time Course data, TRUE IF SS
    public ReadData(File datafile) throws FileNotFoundException, IOException{
        
        try {
            
            BufferedReader input = new BufferedReader (
                                        new InputStreamReader(
                                            new FileInputStream(datafile)));
            
            ArrayList lines = new ArrayList();            
            String line;
            while ((line = input.readLine()) != null) {
                lines.add(line);
            }
            
            String firstline = (String) lines.get(0);
            String[] temp = firstline.split("\\s");
            if(temp[0].equals("Time")){
                SSorTC=false;
            }
            
            if (SSorTC==true){
                IDs= firstline.split("\\s");
                
                String secondline = (String) lines.get(1);
                String[] temp2 = secondline.split("\\s");
                SSvalues = new double[temp2.length];
                
                for (int j =0; j<temp2.length; j++){
                    SSvalues[j]=Double.parseDouble(temp2[j]);
                }
                
                for (int i =0; i<IDs.length;i++){
                    SSdata.put(IDs[i],SSvalues[i]);
                }
            }else{
                IDs=firstline.split("\\s");
                int row = lines.size();
                int col = IDs.length;
                
                while ((line = input.readLine()) != null) {
                    row++;
                }
                
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
                        timeCourse[j-1] = Double.parseDouble(timeCoursedata[j][i]);
                    }                 
                    TCdata.put(IDs[i], timeCourse);
                }                
            }
            
        } catch (Exception e) {
        }
        
    }
    public String[] getID(){
        return IDs;
    }
    
    public double[] getSSvalues(){
        return SSvalues;
    }
    
    public HashMap getSSdata(){
        return SSdata;
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
        
}

