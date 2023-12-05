#include <stdio.h>
#include <math.h>
#include "cgp.h"

#define NUM_INPUTS 2
#define NUM_OUTPUTS 1
#define NUM_SAMPLES 4
#define MAX_INPUTS_PER_GATE 2

double andFunction(const int numInputs, const double *inputs, const double *connectionWeights)
{
    double output = 1;
    for (int i = 0; i < numInputs; i++)
    {
        output *= inputs[i];
    }
    return output;
}
double orFunction(const int numInputs, const double *inputs, const double *connectionWeights)
{
    double output = 0;
    for (int i = 0; i < numInputs; i++)
    {
        output += inputs[i];
    }
    if (output > 1)
        output = 1;
    return output;
}
double notFunction(const int numInputs, const double *inputs, const double *connectionWeights)
{
    return inputs[0] > 0 ? 0 : 1;
}

int main(void)
{

    // PREPARA O DATASET
    int i;

    struct dataSet *data = NULL;

    double inp[NUM_SAMPLES][NUM_INPUTS];
    double out[NUM_SAMPLES][NUM_OUTPUTS];

    inp[0][0] = 0;
    inp[0][1] = 0;
    out[0][0] = 0;

    inp[1][0] = 0;
    inp[1][1] = 1;
    out[1][0] = 1;

    inp[2][0] = 1;
    inp[2][1] = 0;
    out[2][0] = 0;

    inp[3][0] = 1;
    inp[3][1] = 1;
    out[3][0] = 0;

    data = initialiseDataSetFromArrays(NUM_INPUTS, NUM_OUTPUTS, NUM_SAMPLES, inp[0], out[0]);

    saveDataSet(data, "symbolic.data");

    freeDataSet(data);

    ///////////////////
    // TREINA E IMPRIME UMA SOLUÇÃO

    struct parameters *params = NULL;
    struct dataSet *trainingData = NULL;
    struct chromosome *chromo = NULL;

    int numInputs = NUM_INPUTS;
    int numNodes = 15;
    int numOutputs = NUM_OUTPUTS;
    int nodeArity = 2;

    int numGens = 10000;
    int updateFrequency = 500;
    double targetFitness = 0.1;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addCustomNodeFunction(params, andFunction, "and", MAX_INPUTS_PER_GATE);
    addCustomNodeFunction(params, orFunction, "or", MAX_INPUTS_PER_GATE);
    addCustomNodeFunction(params, notFunction, "not", 1);

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    printParameters(params);

    trainingData = initialiseDataSetFromFile("symbolic.data");

    chromo = runCGP(params, trainingData, numGens);

    removeInactiveNodes(chromo);
    printChromosome(chromo, 0);
    saveChromosomeDot(chromo, 0, "chromo.dot");
    saveChromosomeLatex(chromo, 0, "chromo.tex");

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);

    return 0;
}