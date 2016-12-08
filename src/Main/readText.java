/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 *
 * @author chuanfuyap
 */
public class readText {
    
    private double[] fakeparameters;
    
    public readText(File textfile) throws FileNotFoundException, IOException{
        int counter = 0;
        ArrayList<Double> number = new ArrayList<Double>();
        try {
            BufferedReader input = new BufferedReader (
                                        new InputStreamReader(
                                            new FileInputStream(textfile)));

            String line;      //reads lines.     
            
            
            while((line=input.readLine())!=null) {
                counter++;
                number.add(Double.parseDouble(line));
            }
            
            
        } catch (Exception e) {
        }
        fakeparameters=new double[counter];
        for (int i =0; i<number.size(); i++){
            fakeparameters[i]=number.get(i);
        }
    }
    
    public double[] getpara(){
        return fakeparameters;
    }
    
}
