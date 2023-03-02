#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgp.h"

#define PREPARE_DATA
#define RUN_CGP

#define NUM_NODES 300
#define MUTATION_RATE 0.05
#define MAX_GENERATIONS 50000

#define ADDER_SIZE 4
#define NUM_SAMPLES 16 * 16

#define NUM_INPUTS ADDER_SIZE * 2  // 8
#define NUM_OUTPUTS ADDER_SIZE + 1 // 5
#define MAX_INPUTS_PER_GATE 3

double inputs[NUM_SAMPLES][NUM_INPUTS];
double outputs[NUM_SAMPLES][NUM_OUTPUTS];

void generateAdder()
{
    int linha = 0;
    int maxN = pow(2, ADDER_SIZE);
    for (int nA = 0; nA < maxN; nA++)
    {
        for (int nB = 0; nB < maxN; nB++)
        {
            printf("%d >> %d (", linha, nA);
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

double majority(const int numInputs, const double *inputs, const double *connectionWeights)
{
    int nHighs = 0;

    for (int i = 0; i < numInputs; i++)
    {
        if (inputs[i] > 0)
            nHighs++;
        if (nHighs > (numInputs / 2))
            return 1.0;
    }
    return 0.0;
}

double checkTruthTable(struct parameters *params, struct chromosome *chromo, struct dataSet *data)
{
    double linesErrors = 0;
    for (int i = 0; i < getNumDataSetSamples(data); i++)
    {
        int errors = 0;
        executeChromosome(chromo, getDataSetSampleInputs(data, i));

        for (int j = 0; j < getNumChromosomeOutputs(chromo); j++)
        {
            if (getDataSetSampleOutput(data, i, j) != getChromosomeOutput(chromo, j))
                errors += 1;
        }
        if (errors > 0)
            linesErrors += 1;
    }
    return linesErrors;
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

    params = initialiseParameters(NUM_INPUTS, NUM_NODES, NUM_OUTPUTS, MAX_INPUTS_PER_GATE);

    addNodeFunction(params, "and,nand,or,nor,not,xor,xand");
    addCustomNodeFunction(params, majority, "maj", -1);
    setTargetFitness(params, 1);
    setCustomFitnessFunction(params, checkTruthTable, "CTT");

    setEvolutionaryStrategy(params, '+');
    setMu(params, 5);      // The number of parents selected each iteration.
    setLambda(params, 20); // Size of the population.
    // lambda / mu: Number of children generated from each selected parent.
    // (mu, lambda)-ES: A version of evolution strategies where children replace parents.
    // (mu + lambda)-ES: A version of evolution strategies where children and parents are added to the population.
    setMutationRate(params, MUTATION_RATE);
    setUpdateFrequency(params, 500);
    setNumThreads(params, 12);

    printParameters(params);

    trainingData = initialiseDataSetFromFile("symbolic.data");

    chromo = runCGP(params, trainingData, MAX_GENERATIONS);

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