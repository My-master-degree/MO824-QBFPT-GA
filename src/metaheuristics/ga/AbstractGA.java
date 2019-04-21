package metaheuristics.ga;

import java.util.ArrayList;
import java.util.Random;

import problems.Evaluator;
import solutions.Solution;

/**
 * Abstract class for metaheuristic GA (Genetic Algorithms). It consider the
 * maximization of the chromosome fitness.
 *
 * @author ccavellucci, fusberti
 * @param <G> Generic type of the chromosome element (genotype).
 * @param <F> Generic type of the candidate to enter the solution (fenotype).
 */
public abstract class AbstractGA<G extends Number, F> {

    public final static int DEFAULT_CROSSOVER = 1;
    public final static int UNIFORM_CROSSOVER = 2;
    public int CROSSOVER_TYPE;

    @SuppressWarnings("serial")
    public class Population extends ArrayList<Chromosome<F>> {
    }

    /**
     * flag that indicates whether the code should print more information on
     * screen
     */
    public static boolean verbose = true;

    /**
     * a random number generator
     */
    public static final Random rng = new Random(0);

    /**
     * the objective function being optimized
     */
    protected Evaluator<F> ObjFunction;

    /**
     * tempo para execução do GRASP.
     */
    protected Integer tempoExecucao;

    /**
     * Quantidade de geracoes sem melhora para considerar a convergễncia do
     * algoritmo.
     */
    protected Integer geracoesConvengencia;

    /**
     * the size of the population
     */
    protected int popSize;

    /**
     * the size of the chromosome
     */
    protected int chromosomeSize;

    /**
     * the probability of performing a mutation
     */
    protected double mutationRate;

    /**
     * the best solution cost
     */
    protected Double bestCost;

    /**
     * the best solution
     */
    protected Solution<F> bestSol;

    /**
     * the best chromosome, according to its fitness evaluation
     */
    protected Chromosome<F> bestChromosome;

    /**
     * Creates a new solution which is empty, i.e., does not contain any
     * candidate solution element.
     *
     * @return An empty solution.
     */
    public abstract Solution<F> createEmptySol();

    /**
     * A mapping from the genotype (domain) to the fenotype (image). In other
     * words, it takes a chromosome as input and generates a corresponding
     * solution.
     *
     * @param chromosome The genotype being considered for decoding.
     * @return The corresponding fenotype (solution).
     */
    protected abstract Solution<F> decode(Chromosome<F> chromosome);

    protected abstract Chromosome<F> createEmpytChromossome();

    /**
     * Generates a random chromosome according to some probability distribution
     * (usually uniform).
     *
     * @return A random chromosome.
     */
    protected abstract Chromosome<F> generateRandomChromosome();

    /**
     * Mutates a given locus of the chromosome. This method should be preferably
     * called with an expected frequency determined by the
     * {@link #mutationRate}.
     *
     * @param chromosome The genotype being mutated.
     * @param locus The position in the genotype being mutated.
     */
    protected abstract void mutateGene(Chromosome<F> chromosome, Integer locus);

    /**
     * The constructor for the GA class.
     *
     * @param objFunction The objective function being optimized.
     * @param tempoExecucao
     * @param geracoesConvengencia
     * @param popSize Population size.
     * @param mutationRate The mutation rate.
     */
    public AbstractGA(Evaluator<F> objFunction, Integer tempoExecucao, Integer geracoesConvengencia, Integer popSize, Double mutationRate, int crossoverType) {
        this.ObjFunction = objFunction;
        this.tempoExecucao = tempoExecucao;
        this.geracoesConvengencia = geracoesConvengencia;
        this.popSize = popSize;
        this.chromosomeSize = this.ObjFunction.getDomainSize();
        this.mutationRate = mutationRate;
        this.CROSSOVER_TYPE = crossoverType;
    }

    /**
     * The GA mainframe. It starts by initializing a population of chromosomes.
     * It then enters a generational loop, in which each generation goes the
     * following steps: parent selection, crossover, mutation, population update
     * and best solution update.
     *
     * @return The best feasible solution obtained throughout all iterations.
     */
    public Solution<F> solve() {
        long tempoInicial;
        int geracoesSemMelhora = 0;

        /* starts the initial population */
        Population population = initializePopulation();

        bestChromosome = getBestChromosome(population);
        bestSol = decode(bestChromosome);
        System.out.println("(Gen. " + 0 + ") BestSol = " + bestSol);

        /*
         * enters the main loop and repeats until a given number of generations
         */
        tempoInicial = System.currentTimeMillis();
        for (int g = 1; (((System.currentTimeMillis() - tempoInicial) / 1000) / 60) <= this.tempoExecucao && (this.geracoesConvengencia < 1 || geracoesSemMelhora <= this.geracoesConvengencia); g++) {
            geracoesSemMelhora++;

            Population parents = selectParents(population);

            Population offsprings = crossover(parents);

            Population mutants = mutate(offsprings);

            Population newpopulation = selectPopulation(mutants);

            population = newpopulation;

            bestChromosome = getBestChromosome(population);

            if (bestChromosome.getFitnessVal() > bestSol.cost) {
                bestSol = decode(bestChromosome);

                if (verbose) {
                    System.out.println("(Gen. " + g + ") BestSol = " + bestSol);
                }

                geracoesSemMelhora = 0;
            }

        }

        return bestSol;
    }

    /**
     * Randomly generates an initial population to start the GA.
     *
     * @return A population of chromosomes.
     */
    protected Population initializePopulation() {

        Population population = new Population();

        while (population.size() < popSize) {
            Chromosome<F> c = generateRandomChromosome();
            c.calcFitness(ObjFunction);
            population.add(c);
        }

        return population;

    }

    /**
     * Given a population of chromosome, takes the best chromosome according to
     * the fitness evaluation.
     *
     * @param population A population of chromosomes.
     * @return The best chromosome among the population.
     */
    protected Chromosome<F> getBestChromosome(Population population) {

        double bestFitness = Double.NEGATIVE_INFINITY;
        Chromosome<F> bestChromosome = null;
        for (Chromosome<F> c : population) {
            if (c.getFitnessVal() > bestFitness) {
                bestFitness = c.getFitnessVal();
                bestChromosome = c;
            }
        }

        return bestChromosome;
    }

    /**
     * Given a population of chromosome, takes the worst chromosome according to
     * the fitness evaluation.
     *
     * @param population A population of chromosomes.
     * @return The worst chromosome among the population.
     */
    protected Chromosome<F> getWorseChromosome(Population population) {

        double worseFitness = Double.POSITIVE_INFINITY;
        Chromosome<F> worseChromosome = null;
        for (Chromosome<F> c : population) {
            double fitness = c.getFitnessVal();
            if (fitness < worseFitness) {
                worseFitness = fitness;
                worseChromosome = c;
            }
        }

        return worseChromosome;
    }

    /**
     * Selection of parents for crossover using the tournament method. Given a
     * population of chromosomes, randomly takes two chromosomes and compare
     * them by their fitness. The best one is selected as parent. Repeat until
     * the number of selected parents is equal to {@link #popSize}.
     *
     * @param population The current population.
     * @return The selected parents for performing crossover.
     */
    protected Population selectParents(Population population) {

        Population parents = new Population();

        while (parents.size() < popSize) {
            int index1 = rng.nextInt(popSize);
            Chromosome<F> parent1 = population.get(index1);
            int index2 = rng.nextInt(popSize);
            Chromosome<F> parent2 = population.get(index2);
            if (parent1.getFitnessVal() > parent2.getFitnessVal()) {
                parents.add(parent1);
            } else {
                parents.add(parent2);
            }
        }

        return parents;

    }

    protected Population crossover(Population parents) {

        if (this.CROSSOVER_TYPE == UNIFORM_CROSSOVER) {
            return uniformCrossover(parents);
        }

        return defaultCrossover(parents);
    }

    /**
     * The crossover step takes the parents generated by {@link #selectParents}
     * and recombine their genes to generate new chromosomes (offsprings). The
     * method being used is the 2-point crossover, which randomly selects two
     * locus for being the points of exchange (P1 and P2). For example:
     *
     * P1 P2 Parent 1: X1 ... Xi | Xi+1 ... Xj | Xj+1 ... Xn Parent 2: Y1 ... Yi
     * | Yi+1 ... Yj | Yj+1 ... Yn
     *
     * Offspring 1: X1 ... Xi | Yi+1 ... Yj | Xj+1 ... Xn Offspring 2: Y1 ... Yi
     * | Xi+1 ... Xj | Yj+1 ... Yn
     *
     * @param parents The selected parents for crossover.
     * @return The resulting offsprings.
     */
    protected Population defaultCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome<F> parent1 = parents.get(i);
            Chromosome<F> parent2 = parents.get(i + 1);

            int crosspoint1 = rng.nextInt(chromosomeSize + 1);
            int crosspoint2 = crosspoint1 + rng.nextInt((chromosomeSize + 1) - crosspoint1);

            Chromosome<F> offspring1 = createEmpytChromossome();
            Chromosome<F> offspring2 = createEmpytChromossome();

            for (int j = 0; j < chromosomeSize; j++) {
                if (j >= crosspoint1 && j < crosspoint2) {
                    offspring1.add(parent2.get(j));
                    offspring2.add(parent1.get(j));
                } else {
                    offspring1.add(parent1.get(j));
                    offspring2.add(parent2.get(j));
                }
            }

            offspring1.calcFitness(ObjFunction);
            offspring2.calcFitness(ObjFunction);

            offsprings.add(offspring1);
            offsprings.add(offspring2);

        }

        return offsprings;
    }

    protected Population uniformCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {
            Chromosome<F> parent1 = parents.get(i);
            Chromosome<F> parent2 = parents.get(i + 1);

            Chromosome<F> offspring1 = createEmpytChromossome();
            Chromosome<F> offspring2 = createEmpytChromossome();

            for (int j = 0; j < chromosomeSize; j++) {
                if (rng.nextDouble() < 0.5D) {
                    offspring1.add(parent2.get(j));
                    offspring2.add(parent1.get(j));
                } else {
                    offspring1.add(parent1.get(j));
                    offspring2.add(parent2.get(j));
                }
            }

            offspring1.calcFitness(ObjFunction);
            offspring2.calcFitness(ObjFunction);

            offsprings.add(offspring1);
            offsprings.add(offspring2);
        }

        return offsprings;
    }

    /**
     * The mutation step takes the offsprings generated by {@link #crossover}
     * and to each possible locus, perform a mutation with the expected
     * frequency given by {@link #mutationRate}.
     *
     * @param offsprings The offsprings chromosomes generated by the
     * {@link #crossover}.
     * @return The mutated offsprings.
     */
    protected Population mutate(Population offsprings) {

        for (Chromosome<F> c : offsprings) {
            boolean teveMutacao = false;

            for (int locus = 0; locus < chromosomeSize; locus++) {
                if (rng.nextDouble() < mutationRate) {
                    mutateGene(c, locus);
                    teveMutacao = true;
                }
            }

            if (teveMutacao) {
                c.calcFitness(ObjFunction);
            }
        }

        return offsprings;
    }

    /**
     * Updates the population that will be considered for the next GA
     * generation. The method used for updating the population is the elitist,
     * which simply takes the worse chromosome from the offsprings and replace
     * it with the best chromosome from the previous generation.
     *
     * @param offsprings The offsprings generated by {@link #crossover}.
     * @return The updated population for the next generation.
     */
    protected Population selectPopulation(Population offsprings) {

        Chromosome<F> worse = getWorseChromosome(offsprings);
        if (worse.getFitnessVal() < bestChromosome.getFitnessVal()) {
            offsprings.remove(worse);
            offsprings.add(bestChromosome);
        }

        return offsprings;
    }

}
