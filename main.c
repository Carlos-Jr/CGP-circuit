#include <stdio.h>
#include <math.h>
#include "cgp.h"

#define ADDER_SIZE 4
#define NUM_SAMPLES 16

#define NUM_INPUTS ADDER_SIZE * 2  // 8
#define NUM_OUTPUTS ADDER_SIZE + 1 // 5
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

double inputs[NUM_SAMPLES][NUM_INPUTS];
double outputs[NUM_SAMPLES][NUM_OUTPUTS];

void generateAdder()
{
    for (int n = 0; n < NUM_SAMPLES; n++)
    {
        for (int i = 0; i < NUM_INPUTS / 2; i++)
        {
            // Entrada A
            inputs[n][i] = n >> i & 1;
            printf("[%d][%d]:%d", n, i, inputs[n][i]);
            // Entrada B
            inputs[n][i + NUM_INPUTS / 2] = n >> i & 1;
        }
        int calcOut = n * 2;
        for (int i = 0; i < NUM_OUTPUTS; i++)
        {
            outputs[n][i] = calcOut >> i & 1;
        }
    }
}

int main(void)
{
    /////////////////////
    // PREPARA O DATASET
    /////////////////////
    int i;
    struct dataSet *data = NULL;
    generateAdder();
    data = initialiseDataSetFromArrays(NUM_INPUTS, NUM_OUTPUTS, NUM_SAMPLES, inputs[0], outputs[0]);
    saveDataSet(data, "symbolic.data");
    freeDataSet(data);

    /////////////////////
    // RODAR O CGP
    /////////////////////

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