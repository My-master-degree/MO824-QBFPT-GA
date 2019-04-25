package problems.qbf.solvers;

import java.io.IOException;
import metaheuristics.ga.AbstractGA;
import metaheuristics.ga.Chromosome;
import problems.qbf.QBF;
import solutions.Solution;

/**
 * Metaheuristic GA (Genetic Algorithm) for obtaining an optimal solution to a
 * QBF (Quadractive Binary Function -- {@link #QuadracticBinaryFunction}).
 *
 * @author ccavellucci, fusberti
 */
public class GA_QBF extends AbstractGA<Integer, Integer> {

    /**
     * Constructor for the GA_QBF class. The QBF objective function is passed as
     * argument for the superclass constructor.
     *
     * @param tempoExecucao
     * @param geracoesConvengencia
     * @param popSize Size of the population.
     * @param mutationRate The mutation rate.
     * @param filename Name of the file for which the objective function
     * parameters should be read.
     * @throws IOException Necessary for I/O operations.
     */
    public GA_QBF(Integer tempoExecucao, Integer geracoesConvengencia, Integer popSize, Double mutationRate, String filename, int crossoverType, int mutationType) throws IOException {
        super(new QBF(filename), tempoExecucao, geracoesConvengencia, popSize, mutationRate, crossoverType, mutationType);
    }

    /**
     * {@inheritDoc}
     *
     * This createEmptySol instantiates an empty solution and it attributes a
     * zero cost, since it is known that a QBF solution with all variables set
     * to zero has also zero cost.
     */
    @Override
    public Solution<Integer> createEmptySol() {
        Solution<Integer> sol = new Solution<Integer>();
        sol.cost = 0.0;
        return sol;
    }


    /*
	 * (non-Javadoc)
	 * 
	 * @see metaheuristics.ga.AbstractGA#generateRandomChromosome()
     */
    @Override
    protected Chromosome<Integer> generateRandomChromosome() {

        Chromosome<Integer> chromosome = createEmpytChromossome();
        for (int i = 0; i < chromosomeSize; i++) {
            chromosome.add(rng.nextInt(2));
        }

        return chromosome;
    }

    /*
	 * (non-Javadoc)
	 * 
	 * @see
	 * metaheuristics.ga.AbstractGA#mutateGene(metaheuristics.ga.AbstractGA.
	 * Chromosome, java.lang.Integer)
     */
    @Override
    protected void mutateGene(Chromosome<Integer> chromosome, Integer locus) {

        chromosome.set(locus, 1 - chromosome.get(locus));

    }

    /*
	 * (non-Javadoc)
	 * 
	 * @see metaheuristics.ga.AbstractGA#decode(metaheuristics.ga.AbstractGA.
	 * Chromosome)
     */
    @Override
    protected Solution<Integer> decode(Chromosome<Integer> chromosome) {
        Solution<Integer> solution = createEmptySol();

        for (int locus = 0; locus < chromosome.size(); locus++) {
            if (chromosome.get(locus) == 1) {
                solution.add(new Integer(locus));
            }
        }

        ObjFunction.evaluate(solution);
        return solution;
    }

    @Override
    protected Chromosome<Integer> createEmpytChromossome() {
        return new ChromossomeQBF();
    }

    @Override
    protected void endGenerationAction() {
        // TODO Auto-generated method stub

    }

    @Override
    public Boolean mutationCriteria() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     * A main method used for testing the GA metaheuristic.
     *
     */
    public static void main(String[] args) throws IOException {

        long startTime = System.currentTimeMillis();
        GA_QBF ga = new GA_QBF(30, 1000, 100, 1.0 / 100.0, "instances/qbf020", AbstractGA.DEFAULT_CROSSOVER, AbstractGA.DEFAULT_MUTATION);
        Solution<Integer> bestSol = ga.solve();
        System.out.println("maxVal = " + bestSol);
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

    }

}
