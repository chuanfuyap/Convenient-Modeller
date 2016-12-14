#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include </home/chuanfuyap/SOSLib/buildSOSlib/src/sbmlsolver/odeSolver.h>
#include </home/chuanfuyap/SOSLib/buildSOSlib/src/sbmlsolver/sbmlResults.h>
#include </home/chuanfuyap/SOSLib/buildSOSlib/src/sbmlsolver/solverError.h>

typedef struct rowarray{    //struct containing array of metconc
    double* vals;
} rowarray;

void TCrunSolverPK (const char* filelink, const char* variable, int numVals,const double* newparameter,  const char* paranames,  int paramcount, const char* reactionID, int rIDsize, const char* reactionID2,const double* time, int last, rowarray** row, int* rowcount, const char* proteinID, int pIDsize)
{
    int j; 
    int i;
        
    *rowcount= (last+1);
    *row = (rowarray*)malloc(sizeof(rowarray) * (last+1));
    memset(*row, 0, sizeof(rowarray) * (last+1));
    for(i=0;i<(last+1);i++){
        (*row)[i].vals = (double*)malloc(sizeof(double) * (numVals+rIDsize+pIDsize));
        memset((*row)[i].vals, 0, sizeof(double) * (numVals+rIDsize+pIDsize));
    }
    
    //making sure the variables passed isn't empty
    assert(NULL != variable);
    assert(NULL != paranames);
    assert(NULL != reactionID);
    int offset = 0;
    
    char** specieslist =malloc((numVals-1) * sizeof(char*));//freed
    char** paralist = malloc(paramcount * sizeof(char*));//freed
    char** rID = malloc(paramcount * sizeof(char*));//freed
    char** rID2 = malloc(rIDsize * sizeof(char*));//freed
    char** pID = malloc(pIDsize * sizeof(char*));//freed
    
    for(j=0; j<paramcount; j++){
      paralist[j]=malloc(sizeof(char*));
      rID[j]=malloc(sizeof(char*));
    }
    for(j=0; j<rIDsize; j++){      
      rID2[j]=malloc(sizeof(char*));
    }
    for(j=0; j<pIDsize; j++){      
      pID[j]=malloc(sizeof(char*));
    }
    offset = 0;
    for(j=0; j<paramcount; j++){
      paralist[j]=paranames + offset;
      offset += strlen(paralist[j]) + 1;
    }
    offset = 0;
    for(j=0; j<paramcount; j++){
      rID[j]=reactionID + offset;
      offset += strlen(rID[j]) + 1;
    }
    offset = 0;
    for(j=0; j<rIDsize; j++){
      rID2[j]=reactionID2 + offset;
      offset += strlen(rID2[j]) + 1;
    }
    offset = 0;
    for(j=0; j<pIDsize; j++){
      pID[j]=proteinID + offset;
      offset += strlen(pID[j]) + 1;
    }
    
    for(j=0; j<(numVals-1); j++){
      specieslist[j]=malloc(sizeof(char*));
    }
    offset = 0;
    for(j=0; j<(numVals-1); j++){
      specieslist[j]=variable + offset;
      offset += strlen(specieslist[j]) + 1;
    }
    
    SBMLReader_t *sr = SBMLReader_create(); //freed
    SBMLDocument_t *d = SBMLReader_readSBML(sr, filelink); //freed
    Model_t *sbml = SBMLDocument_getModel(d); //freed
    SBMLReader_free(sr); 
    
    //make it so that it has a time step of 1 second, making it easier for me to extract them
    cvodeSettings_t *options = CvodeSettings_createWithTime(time[last], time[last]); //freed        
    
    varySettings_t *parameters = VarySettings_allocate(paramcount,1);
    
    for(i=0; i<paramcount;i++){
        VarySettings_setName(parameters, i, paralist[i], rID[i]);
        VarySettings_setValue(parameters, 0, i, newparameter[i]);
    }
    
    SBMLResultsArray_t *arrayresults = Model_odeSolverBatch(sbml, options, parameters); //freed
    
    SBMLResults_t *results = SBMLResultsArray_getResults(arrayresults, 0);
        
    int tick;
    
    if ( results != NULL )
    {
    for(i=0; i< (last+1); i++){
        //row = time points
        (*row)[i].vals[0]=time[i];
        for ( j=1; j<numVals; j++ ){ //numVals is number of variables+1 (to include time) so since j starts at 1, i dont have to + 1 on numvals
            //column =met conc;
            //specieslist is j-1 to start from 0;
            (*row)[i].vals[j] = TimeCourse_getValue(Species_getTimeCourse(Model_getSpeciesById(sbml, specieslist[(j-1)]), results), time[i]);            
        }//maybe i should consider using species ID to for the timecourse instead of the species object
        tick=0;
        for(j=(numVals);j<(numVals+pIDsize);j++){ //previously it stopped 1 number before numVals so now it starts at numVals and will stop at numVals+rIDsize
            (*row)[i].vals[j] = TimeCourse_getValue(SBMLResults_getTimeCourse(results, pID[tick]), time[i]);
            tick++;
        }
        tick=0;
        for(j=(numVals+pIDsize);j<(numVals+rIDsize+pIDsize);j++){ //previously it stopped 1 number before numVals so now it starts at numVals and will stop at numVals+rIDsize
            (*row)[i].vals[j] = TimeCourse_getValue(SBMLResults_getTimeCourse(results, rID2[tick]), time[i]);
            tick++;
        }
    }
    }
    
    SolverError_dumpAndClearErrors();
    SBMLResultsArray_free(arrayresults);
    CvodeSettings_free(options);
    VarySettings_free(parameters);
    SBMLDocument_free(d);
    
}

void TCarrayCleanupPK(rowarray* row, int count)
{
    int i;
    for( i =0; i < count;i++){
        free(row[i].vals);
    }
    free(row);
}

