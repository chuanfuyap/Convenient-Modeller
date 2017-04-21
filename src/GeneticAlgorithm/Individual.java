/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import java.util.HashMap;
import java.util.Random;
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
    private Random generator = new Random();
    private HashMap metConc;
    private HashMap flux;

    public Individual(int size) {
        this.size = size;
        x = new double[size];
        f = 1e100;
    }
    
    public void initialize(SystemToSolve bioSystem){
        int[] fRIndex=bioSystem.givefRIndex();
        int[] rRIndex=bioSystem.giverRIndex();
        
        int counter=0;
        int counter2=0;
        
        for (int j = 0; j < size; j++) {
            
            if(counter<fRIndex.length && j==fRIndex[counter]) {       //for forward rate constant and backward rate constant range which is thought to have a smaller range in general
                double num = this.RN(-3, 8); // 3 8 are values extracted from SABIO-RK
                x[j]=Math.pow(10, num);
                counter++;
            }else if (counter2<rRIndex.length && j==rRIndex[counter2]) {
                double num=this.RN(-3, 1); // -3 -1 are values extracted from SABIO-RK
                x[j]=Math.pow(10, num);
                counter2++;
            }else{
                double num = this.RN(-5, 5); //
                x[j]=Math.pow(10, num);
            }
        }
    }
    
    private double RN(double lower, double upper) {
        return random.nextUniform(lower, upper);
    }
    
    private double RN() {
        return generator.nextDouble();
    }
    
    public double getF(){
        return f;
    }
    
    public void setF(double fitness){
        f=fitness;
    }
    
}
