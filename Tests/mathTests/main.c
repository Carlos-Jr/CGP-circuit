#include <stdio.h>
#include <math.h>
#include "cgp.h"

#define NUM_INPUTS 1
#define NUM_OUTPUTS 1
#define NUM_SAMPLES 101
#define INPUTRANGE 10.0

double symbolicEq1(double x)
{
    return 2 * x;
}

int main(void)
{

    // PREPARA O DATASET
    int i;

    struct dataSet *data = NULL;

    double inputs[NUM_SAMPLES][NUM_INPUTS];
    double outputs[NUM_SAMPLES][NUM_OUTPUTS];

    double inputTemp;
    double outputTemp;

    for (i = 0; i < NUM_SAMPLES; i++)
    {

        inputTemp = (i * (INPUTRANGE / (NUM_SAMPLES - 1))) - INPUTRANGE / 2;
        outputTemp = symbolicEq1(inputTemp);

        inputs[i][0] = inputTemp;
        outputs[i][0] = outputTemp;
    }

    data = initialiseDataSetFromArrays(NUM_INPUTS, NUM_OUTPUTS, NUM_SAMPLES, inputs[0], outputs[0]);

    saveDataSet(data, "symbolic.data");

    freeDataSet(data);

    ///////////////////
    // TREINA E IMPRIME UMA SOLUÇÃO

    struct parameters *params = NULL;
    struct dataSet *trainingData = NULL;
    struct chromosome *chromo = NULL;

    int numInputs = 1;
    int numNodes = 15;
    int numOutputs = 1;
    int nodeArity = 2;

    int numGens = 10000;
    int updateFrequency = 500;
    double targetFitness = 0.1;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,sub,mul,div,sin");

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