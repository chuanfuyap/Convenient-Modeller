/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SteadyState;

import java.util.ArrayList;

/**
 *
 * @author chuanfuyap
 */
public class CKFunction {
    private double N, D; //numerator and denominator
    private double N2, D2; //numerator and denominator for Partial Derivate
    private int numSub, numProd, numKP, regulation;
    private double Enz, Vf, Vr, KMod, modifier;
    private double[] KSub, KProd;
    private double freg = 1;
    private Double[] Sub, Prod, KiPa;
    private ArrayList substoichio, prodstoichio;
     
    public CKFunction (ArrayList<Double> Substrates, ArrayList<Double> Products, ArrayList<Double> Parameters, double E, int regulation, double modifier, ArrayList substoichio, ArrayList prodstoichio){//E = Enzyme
        this.substoichio = substoichio;
        this.prodstoichio = prodstoichio;
        this.modifier=modifier;
        this.regulation=regulation;
        
        numSub = Substrates.size();
        numProd = Products.size();
        numKP = Parameters.size();
        
        Sub = Substrates.toArray(new Double[numSub]);
        Prod = Products.toArray(new Double[numProd]);
        KiPa = Parameters.toArray(new Double[numKP]);
        Enz = E;
        
        Vf = KiPa[0]; Vr = KiPa[1];
        KSub = new double[numSub];
        KProd = new double[numProd];
        int track = 2;
        if(regulation==2||regulation==3){
            KMod=KiPa[track];
            track++;
        }
        for(int i =0; i<numSub;i++){
            KSub[i]=KiPa[track];
            track++;
        }
        for(int i =0; i<numProd;i++){
            KProd[i]=KiPa[track];
            track++;
        }
        
        if(regulation == 2){
            freg= 1 + (modifier/KMod);
        }
        
        if(regulation == 3){
            freg= KMod / (KMod+ modifier);
        }
        
        double topsub = 1;
        for (int i =0; i <numSub; i ++){
            topsub*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
        }
        double topprod = 1;
        for(int i =0; i<numProd; i++){
            topprod*=Math.pow((Prod[i]/KProd[i]), (int)prodstoichio.get(i));
        }

        N = E * freg * ((Vf*topsub)-(Vr*topprod));
        
        
        double btmsub=1;
        for(int i =0; i<numSub;i++){
            double sub=0;
            for(int j=0; j<(int)substoichio.get(i)+1;j++){
                sub+=Math.pow((Sub[i]/KSub[i]),j);
            }
            btmsub*=sub;
        }
        
        double btmprod=1;
        for(int i =0; i<numProd;i++){
            double prod=0;
            for(int j=0;j<(int)prodstoichio.get(i)+1;j++){
                prod+=Math.pow((Prod[i]/KProd[i]),j);

            }
            btmprod*=prod;
        }
        
        D = btmsub + btmprod - 1;
    }
    
    public double SolveFunction(){
//        System.out.println(N+"\t"+D);
        return (N/D);
    }
    
    public double getPartialDerivate(boolean subORprod, int index){ //boolean subORDprod, true if its substrate, false if product
        double dN, dD;
        
        if(subORprod==true){
            dN=Enz*freg*Vf;
            dD=1;
            for(int i=0; i < numSub;i++){
                if(i==index){
                    double temp1 = Math.pow(Sub[i], (int)substoichio.get(i)-1);
                    double temp2 = Math.pow(KSub[i], (int)substoichio.get(i));
                    dN*=((int)substoichio.get(i))*(temp1/temp2);
                    
                    double sub=0;
                    for(int j=1; j<(int)substoichio.get(i)+1;j++){
                        sub+=j*(Math.pow((Sub[i]/KSub[i]),j-1));
                    }
                    dD*=sub;
                }else{
                    dN*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
                    
                    double sub=0;
                    for(int j=0; j<(int)substoichio.get(i)+1;j++){              //because these becomes a constant in the differentiation
                        sub+=Math.pow((Sub[i]/KSub[i]),j);                      //so i can just add them up and multiply to whatever its differentiated
                    }
                    dD*=sub;
                }
            }
            N2 = (D*dN) - (N*dD);
            D2 = Math.pow(D, 2);
        }else{
            dN=Enz*freg*Vr;
            dD=1;
            for(int i=0; i < numProd;i++){
                if(i==index){
                    double temp1 = Math.pow(Prod[i], (int)prodstoichio.get(i)-1);
                    double temp2 = Math.pow(KProd[i], (int)prodstoichio.get(i));
                    dN*=-((int)prodstoichio.get(i))*(temp1/temp2);
                    
                    double prod=0;
                    for(int j=1; j<(int)prodstoichio.get(i)+1;j++){
                        prod+=j*(Math.pow((Prod[i]/KProd[i]),j-1));
                    }
                    dD*=prod;
                }else{
                    dN*=Math.pow((Prod[i]/KProd[i]),(int)prodstoichio.get(i));
                    
                    double prod=0;
                    for(int j=0; j<(int)prodstoichio.get(i)+1;j++){             //because these becomes a constant in the differentiation
                        prod+=Math.pow((Prod[i]/KProd[i]),j);                   //so i can just add them up and multiply to whatever its differentiated
                    }
                    dD*=prod;
                }
            }            
            N2 = (D*dN) - (N*dD);
            D2 = Math.pow(D, 2);
        }
        return N2/D2;
    }
    
    public double getPartialDerivateWithMod(){
        double value=0;
        double newN=1;
        double newD=D;
        
        double topsub = 1;
        for (int i =0; i <numSub; i ++){
            topsub*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
        }
        double topprod = 1;
        for(int i =0; i<numProd; i++){
            topprod*=Math.pow((Prod[i]/KProd[i]), (int)prodstoichio.get(i));
        }
        newN*= ((Vf*topsub)-(Vr*topprod));
        
        if(regulation==2){
            value=(Enz*newN)/(KMod*newD);
        }
        if(regulation==3){
            double top=-newD*KMod*Enz*newN;
            double btm=(KMod*newD)+(modifier*newD);
            double btm2=Math.pow(btm,2);
            
            value=top/btm;
            
        }
        
        return value;
    }

    public double getPartialDerivateWithModSubProd(boolean subORprod, int index) {
        double newN=1;
        double newD=D;
        
        double dN, dD, ddN;
        
        if(regulation==2){
                if(subORprod==true){
                    dN=Enz;
                    ddN=Vf;

                    dD=1;
                    for(int i =0; i<numSub;i++){
                        if(i==index){
                            double temp1 = Math.pow(Sub[i], (int)substoichio.get(i)-1);
                            double temp2 = Math.pow(KSub[i], (int)substoichio.get(i));
                            ddN*= ((int)substoichio.get(i))*(temp1/temp2);                    
                        }else{
                            ddN*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
                        }
                    }

                    for(int i=0; i < numSub;i++){

                        if(i==index){
                            double temp1 = Math.pow(Sub[i], (int)substoichio.get(i));
                            double temp2 = Math.pow(KSub[i], (int)substoichio.get(i));
                            dN*=((int)substoichio.get(i)+1)*(temp1/temp2);
                            dN*=(Vf/KMod);
                            dN+=ddN;

                            double sub=0;
                            for(int j=1; j<(int)substoichio.get(i)+1;j++){
                                sub+=j*(Math.pow((Sub[i]/KSub[i]),j-1));
                            }
                            dD*=sub;
                        }else{
                            dN*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));                    

                            double sub=0;
                            for(int j=0; j<(int)substoichio.get(i)+1;j++){
                                sub+=Math.pow((Sub[i]/KSub[i]),j);
                            }
                            dD*=sub;
                        }
                    }
                    double nuprod=Vr/KMod;
                    for(int i=0; i < numProd;i++){
                        nuprod*=Math.pow((Prod[i]/KProd[i]),(int)prodstoichio.get(i));
                    }
                    dN-=nuprod;
                    N2 = (D*dN) - (N*dD);
                    D2 = Math.pow(D, 2);
                }
                else{       //product;
                    dN=Enz;
                    ddN=Vr;

                    dD=1;
                    for(int i =0; i<numProd;i++){
                        if(i==index){
                            double temp1 = Math.pow(Prod[i], (int)prodstoichio.get(i)-1);
                            double temp2 = Math.pow(KProd[i], (int)prodstoichio.get(i));
                            ddN*= ((int)prodstoichio.get(i))*(temp1/temp2);                 
                        }else{
                            ddN*=Math.pow((Prod[i]/KProd[i]),(int)prodstoichio.get(i));
                        }
                    }

                    for(int i=0; i < numProd;i++){
                        if(i==index){
                            double temp1 = Math.pow(Prod[i], (int)prodstoichio.get(i));
                            double temp2 = Math.pow(KProd[i], (int)prodstoichio.get(i));
                            dN*=-((int)prodstoichio.get(i)+1)*(temp1/temp2);
                            dN*=(Vr/KMod);
                            dN-=ddN;

                            double prod=0;
                            for(int j=1; j<(int)prodstoichio.get(i)+1;j++){
                                prod+=j*(Math.pow((Prod[i]/KProd[i]),j-1));
                            }
                            dD*=prod;
                        }else{
                            dN*=Math.pow((Prod[i]/KProd[i]),(int)prodstoichio.get(i));

                            double prod=0;
                            for(int j=0; j<(int)prodstoichio.get(i)+1;j++){
                                prod+=Math.pow((Prod[i]/KProd[i]),j);
                            }
                            dD*=prod;
                        }
                    }
                    double nusub=Vf/KMod;
                    for(int i =0; i<numSub; i++){
                        nusub*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
                    }
                    dN+=nusub;
                    N2 = (D*dN) - (N*dD);
                    D2 = Math.pow(D, 2);
                }
        
        }
        
        else if (regulation==3){      //inhibitor;
            
            newN=N;
            newD=D;
            
            if(subORprod==true){            //inhibitor modifier and a substrate
                dN=1;
                dD=1;
                
                //first solve differentiation on the numerator dN;
                double uN = Enz*KMod;
                double vN = KMod+modifier;
                double duN = uN*Vf;
                double dvN = 1;
                for(int i=0; i < numSub;i++){
                    if(i==index){
                        double temp1 = Math.pow(Sub[i], (int)substoichio.get(i)-1);
                        double temp2 = Math.pow(KSub[i], (int)substoichio.get(i));
                        duN*=((int)substoichio.get(i))*(temp1/temp2);                        
                    }else{
                        duN*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));                   
                    }
                }                
                double nusub=Vf;
                for(int i=0; i < numSub;i++){
                    nusub*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
                }
                double nuprod=Vr;
                for(int i=0; i < numProd;i++){
                    nuprod*=Math.pow((Prod[i]/KProd[i]), (int)prodstoichio.get(i));
                }
                uN*=(nusub-nuprod);
                double vN2 = Math.pow(vN, 2);
                dN*=((duN*vN)-(dvN*uN))/vN2;
                //done with dN
                
                //now solve differentiation for dD;
                for(int i=0; i < numSub;i++){
                    if(i==index){
                        double sub=0;
                        for(int j=1; j<(int)substoichio.get(i)+1;j++){
                            sub+=j*(Math.pow((Sub[i]/KSub[i]),j-1));
                        }
                        dD*=sub;
                    }else{
                        double sub=0;
                        for(int j=0; j<(int)substoichio.get(i)+1;j++){
                            sub+=Math.pow((Sub[i]/KSub[i]),j);
                        }
                        dD*=sub;
                    }
                }
                N2 = (newD*dN) - (newN*dD);
                D2 = Math.pow(newD, 2);
            }
            
            else{                  //inhibitor modifier and a product
                
                dN=1;
                dD=1;
                
                //first solve differentiation on the numerator dN;
                double uN = Enz*KMod;
                double vN = KMod+modifier;
                double duN = uN*Vr;
                double dvN = 1;
                for(int i=0; i < numProd;i++){
                    if(i==index){
                        double temp1 = Math.pow(Prod[i], (int)prodstoichio.get(i)-1);
                        double temp2 = Math.pow(KProd[i], (int)prodstoichio.get(i));
                        duN*=-((int)prodstoichio.get(i))*(temp1/temp2);                        
                    }else{
                        duN*=Math.pow((Prod[i]/KProd[i]),(int)prodstoichio.get(i));
                    }
                }                
                double nusub=Vf;
                for(int i=0; i < numProd;i++){
                    nusub*=Math.pow((Sub[i]/KSub[i]),(int)substoichio.get(i));
                }
                double nuprod=Vr;
                for(int i=0; i < numSub;i++){
                    nuprod*=Math.pow((Prod[i]/KProd[i]), (int)prodstoichio.get(i));
                }
                uN*=(nusub-nuprod);
                double vN2 = Math.pow(vN, 2);
                double top = (duN*vN)-(dvN*uN);
                dN*=top/vN2;                
                //done with dN
                
                //now solve differentiation for dD;
                for(int i=0; i < numProd;i++){
                    if(i==index){
                        double prod=0;
                        for(int j=1; j<(int)prodstoichio.get(i)+1;j++){
                            prod+=j*(Math.pow((Prod[i]/KProd[i]),j-1));
                        }
                        dD*=prod;
                    }else{
                        double prod=0;
                        for(int j=0; j<(int)prodstoichio.get(i)+1;j++){
                            prod+=Math.pow((Prod[i]/KProd[i]),j);
                        }
                        dD*=prod;
                    }
                }
                N2 = (newD*dN) - (newN*dD);
                D2 = Math.pow(newD, 2);
            }
            
        }
        return N2/D2;
    }
    
}
