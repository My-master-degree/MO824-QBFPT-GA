/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package problems.qbfpt.solvers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import metaheuristics.ga.AbstractGA;
import metaheuristics.ga.Chromosome;
import problems.qbf.solvers.GA_QBF;
import problems.qbfpt.Triple;
import problems.qbfpt.TripleElement;
import solutions.Solution;

/**
 *
 * @author Cristina Bazzano
 * @author Jônatas Trabuco Belotti [jonatas.t.belotti@hotmail.com]
 * @author Matheus Diógenes Andrade
 */
public class GA_QBFPT extends GA_QBF {

    /**
     * List of element objects used in prohibited triples. These objects
     * represents the variables of the model.
     */
    private TripleElement[] tripleElements;

    /**
     * List of prohibited triples.
     */
    private Triple[] triples;

    public GA_QBFPT(Integer tempoExecucao, Integer geracoesConvengencia, Integer popSize, Double mutationRate, String filename, int crossoverType) throws IOException {
        super(tempoExecucao, geracoesConvengencia, popSize, mutationRate, filename, crossoverType);

        generateTripleElements();
        generateTriples();
    }

    @Override
    protected Chromosome<Integer> generateRandomChromosome() {
    	//makeCL();
        Chromosome<Integer> chromosome = createEmpytChromossome();
        
        for (int i = 0; i < chromosomeSize; i++) {
            if (updateCL(chromosome).contains(i)) {
                chromosome.add(rng.nextInt(2));
            } else {
                chromosome.add(0);
            }
        }
        
        System.out.println("chromosome" + chromosome);

        return chromosome;
    }

    @Override
    protected Population defaultCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome<Integer> parent1 = parents.get(i);
            Chromosome<Integer> parent2 = parents.get(i + 1);

//            System.out.println("parent1 "+parent1);
//            System.out.println("parent2 "+parent2);
            
            // encontrar indice a partir do inicio em que parent1 fica diferente de parent2
//            int crossbegin, crossend;
//            for (crossbegin = 0; crossbegin < chromosomeSize; crossbegin++) {
//            	if (parent1.get(crossbegin) != parent2.get(crossbegin))
//            		break;
//            }
//            if (crossbegin == chromosomeSize) // parents are equal
//            {
//            	System.out.println("Parents are equall!!");
//            	continue;
//            }
//            for (crossend = chromosomeSize-1; crossend > crossbegin; crossend--) {
//            	if (parent1.get(crossend) != parent2.get(crossend))
//            		break;
//            }
//            if (crossend == crossbegin)
//            	crossend++;
//            
//            int crosspoint1 = crossbegin + rng.nextInt(crossend + 1 - crossbegin);
//            int crosspoint2 = crosspoint1 + rng.nextInt((crossend + 1) - crosspoint1);
            int crosspoint1 = rng.nextInt(chromosomeSize + 1);
            int crosspoint2 = crosspoint1 + rng.nextInt((chromosomeSize + 1) - crosspoint1);

            Chromosome<Integer> offspring1 = createEmpytChromossome();
            Chromosome<Integer> offspring2 = createEmpytChromossome();

            ArrayList<Integer> CL1 = makeCL();
            ArrayList<Integer> CL2 = makeCL();

            for (int j = 0; j < chromosomeSize; j++) {
                int cand1, cand2;

                if (j >= crosspoint1 && j < crosspoint2) {
                    cand1 = parent2.get(j);
                    cand2 = parent1.get(j);
                } else {
                    cand1 = parent1.get(j);
                    cand2 = parent2.get(j);
                }

                if (!CL1.contains(j)) {
                    cand1 = 0;
                }
                if (!CL2.contains(j)) {
                    cand2 = 0;
                }

                offspring1.add(cand1);
                offspring2.add(cand2);

                CL1 = updateCL(offspring1);
                CL2 = updateCL(offspring2);
            }

            offspring1.calcFitness(ObjFunction);
            offspring2.calcFitness(ObjFunction);

            offsprings.add(offspring1);
            offsprings.add(offspring2);

        }

        return offsprings;
    }

    @Override
    protected Population uniformCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome<Integer> parent1 = parents.get(i);
            Chromosome<Integer> parent2 = parents.get(i + 1);

            Chromosome<Integer> offspring1 = createEmpytChromossome();
            Chromosome<Integer> offspring2 = createEmpytChromossome();

            ArrayList<Integer> CL1 = makeCL();
            ArrayList<Integer> CL2 = makeCL();

            for (int j = 0; j < chromosomeSize; j++) {
                int cand1, cand2;

                if (rng.nextDouble() >= 0.5D) {
                    cand1 = parent2.get(j);
                    cand2 = parent1.get(j);
                } else {
                    cand1 = parent1.get(j);
                    cand2 = parent2.get(j);
                }

                if (!CL1.contains(j)) {
                    cand1 = 0;
                }
                if (!CL2.contains(j)) {
                    cand2 = 0;
                }

                offspring1.add(cand1);
                offspring2.add(cand2);

                CL1 = updateCL(offspring1);
                CL2 = updateCL(offspring2);
            }

            offspring1.calcFitness(ObjFunction);
            offspring2.calcFitness(ObjFunction);

            offsprings.add(offspring1);
            offsprings.add(offspring2);

        }

        return offsprings;
    }

    @Override
    protected void mutateGene(Chromosome<Integer> chromosome, Integer locus) {
        if (chromosome.get(locus) == 1) {
            chromosome.set(locus, 0);
        } else if (updateCL(chromosome).contains(locus)) {
            chromosome.set(locus, 1);
        }
    }

    /**
     * Linear congruent function l used to generate pseudo-random numbers.
     */
    public int l(int pi1, int pi2, int u, int n) {
        return 1 + ((pi1 * u + pi2) % n);
    }

    /**
     * Function g used to generate pseudo-random numbers
     */
    public int g(int u, int n) {
        int pi1 = 131;
        int pi2 = 1031;
        int lU = l(pi1, pi2, u, n);

        if (lU != u) {
            return lU;
        } else {
            return 1 + (lU % n);
        }
    }

    /**
     * Function h used to generate pseudo-random numbers
     */
    public int h(int u, int n) {
        int pi1 = 193;
        int pi2 = 1093;
        int lU = l(pi1, pi2, u, n);
        int gU = g(u, n);

        if (lU != u && lU != gU) {
            return lU;
        } else if ((1 + (lU % n)) != u && (1 + (lU % n)) != gU) {
            return 1 + (lU % n);
        } else {
            return 1 + ((lU + 1) % n);
        }
    }

    /**
     * That method generates a list of objects (Triple Elements) that represents
     * each binary variable that could be inserted into a prohibited triple
     */
    private void generateTripleElements() {
        int n = ObjFunction.getDomainSize();
        this.tripleElements = new TripleElement[n];

        for (int i = 0; i < n; i++) {
            tripleElements[i] = new TripleElement(i);
        }
    }

    /**
     * Method that generates a list of n prohibited triples using l g and h
     * functions
     */
    private void generateTriples() {
        int n = ObjFunction.getDomainSize();
        this.triples = new Triple[ObjFunction.getDomainSize()];

        for (int u = 1; u <= n; u++) {
            TripleElement te1, te2, te3;
            Triple novaTripla;

            te1 = tripleElements[u - 1];
            te2 = tripleElements[g(u - 1, n) - 1];
            te3 = tripleElements[h(u - 1, n) - 1];
            novaTripla = new Triple(te1, te2, te3);

            Collections.sort(novaTripla.getElements(), new Comparator<TripleElement>() {
                public int compare(TripleElement te1, TripleElement te2) {
                    return te1.getIndex().compareTo(te2.getIndex());
                }
            });
            //novaTripla.printTriple();
            this.triples[u - 1] = novaTripla;
        }
    }

    /**
     * A GRASP CL generator for MAXQBFPT problem
     *
     * @return A list of candidates to partial solution
     */
    public ArrayList<Integer> makeCL() {
        int n = ObjFunction.getDomainSize();
        ArrayList<Integer> _CL = new ArrayList<Integer>(n);

        for (TripleElement tripElem : this.tripleElements) {
            tripElem.setAvailable(true);
            tripElem.setSelected(false);
            _CL.add(tripElem.getIndex());
        }

        return _CL;
    }

    /**
     * The CL updater for MAXQBFPT problem
     *
     */
    public ArrayList<Integer> updateCL(Chromosome<Integer> cro) {
        ArrayList<Integer> _CL = new ArrayList<Integer>();
        Solution<Integer> solAtual = decode(cro);

        if (solAtual != null) {
//            for (Integer e : solAtual) {
//                this.tripleElements[e].setSelected(true);
//                this.tripleElements[e].setAvailable(false);
//            }
        	for (int i = 0; i < chromosomeSize; i ++) {
        		if (solAtual.contains(i)) {
        			this.tripleElements[i].setSelected(true);
        			this.tripleElements[i].setAvailable(false);
        		} else {
        			this.tripleElements[i].setSelected(false);
        			this.tripleElements[i].setAvailable(true);
        		}
        	}
        }

        for (Triple trip : this.triples) {
            TripleElement te0, te1, te2;
            te0 = trip.getElements().get(0);
            te1 = trip.getElements().get(1);
            te2 = trip.getElements().get(2);

            if (te0.getSelected() && te1.getSelected()) {
                te2.setAvailable(false);
            } else if (te0.getSelected() && te2.getSelected()) {
                te1.setAvailable(false);
            } else if (te1.getSelected() && te2.getSelected()) {
                te0.setAvailable(false);
            }
        }

        for (TripleElement tripElem : this.tripleElements) {
            if (!tripElem.getSelected() && tripElem.getAvailable()) {
                _CL.add(tripElem.getIndex());
            }
        }

        return _CL;
    }

    /**
     * A main method used for testing the GA metaheuristic.
     *
     */
    public static void main(String[] args) throws IOException {

        long startTime = System.currentTimeMillis();
        GA_QBFPT ga = new GA_QBFPT(30, 1000, 100, 1.0 / 100.0, "instances/qbf040", AbstractGA.DEFAULT_CROSSOVER);
        Solution<Integer> bestSol = ga.solve();
        System.out.println("maxVal = " + bestSol);
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

    }

}
