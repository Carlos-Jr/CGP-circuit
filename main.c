#include <stdio.h>
#include <math.h>
#include "cgp.h"

// #define PREPARE_DATA
#define RUN_CGP

#define ADDER_SIZE 3
#define NUM_SAMPLES 40320

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
    int linha = 0;
    int maxN = pow(2, ADDER_SIZE);
    for (int nA = 0; nA < maxN; nA++)
    {
        for (int nB = nA; nB < maxN; nB++)
        {
            printf("%d (", nA);
            // Entrada A
            for (int i = 0; i < NUM_INPUTS / 2; i++)
            {
                int bitInA = (nA >> ((NUM_INPUTS / 2) - i - 1)) & 1;
                inputs[linha][i] = (double)bitInA;
                printf("%d", (int)inputs[linha][i]);
            }
            printf(") + %d (", nB);
            // Entrada B
            for (int i = 0; i < NUM_INPUTS / 2; i++)
            {
                int bitInB = (nB >> ((NUM_INPUTS / 2) - i - 1)) & 1;
                inputs[linha][i + NUM_INPUTS / 2] = (double)bitInB;
                printf("%d", (int)inputs[linha][i + NUM_INPUTS / 2]);
            }
            int calcOut = nA + nB;
            printf(") = %d (", calcOut);
            for (int i = 0; i < NUM_OUTPUTS; i++)
            {
                int bitInOut = (calcOut >> NUM_OUTPUTS - i - 1) & 1;
                outputs[linha][i] = (double)bitInOut;
                printf("%d", (int)outputs[linha][i]);
            }
            printf(")\n");
            linha++;
        }
    }
}

int main(void)
{
    /////////////////////
    // PREPARA O DATASET
    /////////////////////
#ifdef PREPARE_DATA
    int i;
    struct dataSet *data = NULL;
    generateAdder();
    data = initialiseDataSetFromArrays(NUM_INPUTS, NUM_OUTPUTS, NUM_SAMPLES, inputs[0], outputs[0]);
    saveDataSet(data, "symbolic.data");
    freeDataSet(data);
#endif
    /////////////////////
    // RODAR O CGP
    /////////////////////
#ifdef RUN_CGP

    struct parameters *params = NULL;
    struct dataSet *trainingData = NULL;
    struct chromosome *chromo = NULL;

    int numInputs = NUM_INPUTS;
    int numNodes = 15;
    int numOutputs = NUM_OUTPUTS;
    int nodeArity = 2;

    int numGens = 10000;
    int updateFrequency = 500;
    double targetFitness = 1;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addCustomNodeFunction(params, andFunction, "and", MAX_INPUTS_PER_GATE);
    addCustomNodeFunction(params, orFunction, "or", MAX_INPUTS_PER_GATE);
    addCustomNodeFunction(params, notFunction, "not", 1);

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    setNumThreads(params, 12);

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
#endif
    return 0;
}