/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import org.apache.commons.math3.random.RandomDataGenerator;

/**
 *
 * @author chuanfuyap
 */
public class Individual {

    public int size;            // Number of parameters
    public double[] x;          // The array holding all the parameters for the individual
    private double f;            // Fitness value for the individual;
    private RandomDataGenerator random = new RandomDataGenerator();     //apache maths random number generator for uniform distribution
    private double[] output;

    public Individual(int size) {
        this.size = size;
        x = new double[size];
        f = 1e100;
    }
   
    public double getF(){
        return f;
    }
    
    public void setF(double fitness){
        f=fitness;
    }
    
    public void setOutput(double[] output){
        this.output =output;
    }
    
    public double[] getOutput(){
        return output;
    }
}
