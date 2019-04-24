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
import metaheuristics.ga.AbstractGA.Population;
import problems.qbf.QBF;
import problems.qbf.solvers.ChromossomeQBF;
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
public class GA_QBFPT extends AbstractGA<Integer, Integer> {

    /**
     * List of element objects used in prohibited triples. These objects
     * represents the variables of the model.
     */
    private TripleElement[] tripleElements;

    /**
     * List of prohibited triples.
     */
    private Triple[] triples;
    
    private Double esperanca;
    
    public final static int XOR_CROSSOVER = 3;

    public GA_QBFPT(Integer tempoExecucao, Integer geracoesConvengencia, Integer popSize, Double mutationRate, String filename, int crossoverType, int mutationType) throws IOException {
//        super(tempoExecucao, geracoesConvengencia, popSize, mutationRate, filename, crossoverType);
        super(new QBF(filename), tempoExecucao, geracoesConvengencia, popSize, mutationRate, crossoverType, mutationType);
        esperanca = 0d;
        generateTripleElements();
        generateTriples();
    }

    protected void resetTripleElementsQttUsed() {
    	for (int i = 0; i < tripleElements.length; i++) {
    		tripleElements[i].qttUsed = 0;
		}
    	this.esperanca = 0d;
    }
    
    protected void incrementEsperanca() {
    	esperanca += 1/this.tripleElements.length; 
    }
    
    protected void decrementEsperanca() {
    	esperanca -= 1/this.tripleElements.length; 
    }
    
    @Override
    protected Population crossover(Population parents) {

        if (this.CROSSOVER_TYPE == UNIFORM_CROSSOVER) {
            return uniformCrossover(parents);
        } else if (CROSSOVER_TYPE == XOR_CROSSOVER) {
        	return uniformXorCrossover(parents);
        }

        return defaultCrossover(parents);
    }
    
    @Override
    protected Chromosome<Integer> generateRandomChromosome() {
        Chromosome<Integer> chromosome = createEmpytChromossome();
        ArrayList<Integer> listIndices = new ArrayList<Integer>();

        for (int i = 0; i < chromosomeSize; i++) {
            chromosome.add(0);
            listIndices.add(i);
        }

        Collections.shuffle(listIndices);

        for (int i : listIndices) {
            if (updateCL(chromosome).contains(i)) {
            	Integer used = rng.nextInt(2);
            	if (used == 1) {
            		this.tripleElements[i].qttUsed++;    
            		this.incrementEsperanca();
            	}
                chromosome.set(i, used);
            }
        }

        //System.out.println("chromosome" + chromosome);
        return chromosome;
    }

    protected Integer diffChromosome(Chromosome<Integer> c1, Chromosome<Integer> c2) {
        int diff = 0;

        for (int i = 0; i < chromosomeSize; i++) {
            if (c1.get(i) != c2.get(i)) {
                diff++;
            }
        }

        return diff;
    }

    @Override
    protected Population selectParents(Population population) {

        Population parents = new Population();
        Chromosome<Integer> parent1, parent2;

        while (parents.size() < popSize) {
            // find a parent that is at least 10% different from the last added one
            do {
                int index1 = rng.nextInt(popSize);
                parent1 = population.get(index1);
            } while (parents.size() > 0
                    && diffChromosome(parents.get(parents.size() - 1), parent1) < 0.1 * chromosomeSize);

            // try to find a different parent from at least one element
            do {
                int index2 = rng.nextInt(popSize);
                parent2 = population.get(index2);
            } while (parents.size() > 0
                    && diffChromosome(parents.get(parents.size() - 1), parent2) == 0);

            if (parent1.getFitnessVal() > parent2.getFitnessVal()) {
                parents.add(parent1);
            } else {
                parents.add(parent2);
            }
        }

        return parents;

    }
    
    // return the min dist of the chromossome to any individual of the given population
    // if a clone is found return 0
    protected Integer offspringMinDist(Population offsprings, Chromosome<Integer> chrom) {
    	int mindist = chromosomeSize, dist;
    	
    	for (Chromosome<Integer> individuo : offsprings) {
    		dist = diffChromosome(individuo, chrom);
    		if (dist == 0)
    			return dist; // return when the first clone is found
    		if (dist < mindist)
    			mindist = dist;
    	}
    	
    	return mindist;
    }

    protected ArrayList<Integer> xorPos(Chromosome<Integer> c1, Chromosome<Integer> c2, boolean ones) {
    	ArrayList<Integer> pos = new ArrayList<Integer>();
    	
    	for (int i = 0; i < chromosomeSize; i++) {
    		boolean e1 = (c1.get(i)==1);
    		boolean e2 = (c2.get(i)==1);
    		
    		if (e1 ^ e2)
    		{
    			if (ones)
    				pos.add(i);
    		} else {
    			if (!ones)
    				pos.add(i);
    		}
    	}
    	
    	return pos;    	
    }
    
    @Override
    protected Population defaultCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

            Chromosome<Integer> parent1 = parents.get(i);
            Chromosome<Integer> parent2 = parents.get(i + 1);

            // encontrar indice a partir do inicio em que parent1 fica diferente de parent2
            int crossbegin, crossend;
            for (crossbegin = 0; crossbegin < chromosomeSize; crossbegin++) {
                if (parent1.get(crossbegin) != parent2.get(crossbegin)) {
                    break;
                }
            }
            if (crossbegin == chromosomeSize) // parents are equal
            {
                // mutate parents at random and add them
            	mutateGeneCL(updateCL(parent1), parent1, rng.nextInt(chromosomeSize));

                offsprings.add(parent1);
                offsprings.add(parent2);
                //System.out.println("Parents are equall!!");
                //System.out.println("parent1 "+parent1);
                //System.out.println("parent2 "+parent2);
                continue;
            }
            // encontrar indice a partir do final em que parents ficam diferentes
            for (crossend = chromosomeSize - 1; crossend > crossbegin; crossend--) {
                if (parent1.get(crossend) != parent2.get(crossend)) {
                    break;
                }
            }

            int crosspoint1 = crossbegin + rng.nextInt(crossend + 1 - crossbegin);
            int crosspoint2 = crosspoint1 + rng.nextInt((crossend + 1) - crosspoint1);

//            int crosspoint1 = rng.nextInt(chromosomeSize + 1);
//            int crosspoint2 = crosspoint1 + rng.nextInt((chromosomeSize + 1) - crosspoint1);
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

                if (cand1 == 1 && !CL1.contains(j)) {
                    cand1 = 0;
                }
                if (cand2 == 1 && !CL2.contains(j)) {
                    cand2 = 0;
                }
                
                if (cand1 == 1) {
                	this.tripleElements[j].qttUsed++;
                	this.incrementEsperanca();
                }
                if (cand2 == 1) {
                	this.tripleElements[j].qttUsed++;
                	this.incrementEsperanca();
                }

                offspring1.add(cand1);
                offspring2.add(cand2);

                CL1 = updateCL(offspring1);
                CL2 = updateCL(offspring2);
            }

            offspring1.calcFitness(ObjFunction);
            offspring2.calcFitness(ObjFunction);

            if (offspringMinDist(offsprings, offspring1) == 0) {
            	// existe um individuo igual a offspring1, mutate offspring and add
            	mutateGeneCL(CL1, offspring1, rng.nextInt(chromosomeSize));
            	offspring1.calcFitness(ObjFunction);
            	offsprings.add(offspring1);
            } else {
            	offsprings.add(offspring1);
            }
            if (offspringMinDist(offsprings, offspring2) == 0) {
            	// existe um individuo igual a offspring2, mutate offspring and add
            	mutateGeneCL(CL2, offspring2, rng.nextInt(chromosomeSize));
            	offspring2.calcFitness(ObjFunction);
            	offsprings.add(offspring2);
            } else {
            	offsprings.add(offspring2);
            }            
        }

        return offsprings;
    }


    protected Population uniformXorCrossover(Population parents) {
        Population offsprings = new Population();

        for (int i = 0; i < popSize; i = i + 2) {

//            Chromosome<Integer> parent1 = parents.get(i);
//            Chromosome<Integer> parent2 = parents.get(i + 1);

            Chromosome<Integer> offspring1 = parents.get(i);
            Chromosome<Integer> offspring2 = parents.get(i + 1);

            ArrayList<Integer> CL1 = updateCL(offspring1);
            ArrayList<Integer> CL2 = updateCL(offspring2);
            
            ArrayList<Integer> crossPos = xorPos(offspring1, offspring2, true);
            
            if (crossPos.size() == 0) // parents are equal
            {
                // mutate parents at random and add them
            	mutateGeneCL(CL1, offspring1, rng.nextInt(chromosomeSize));

                offsprings.add(offspring1);
                offsprings.add(offspring2);
                //System.out.println("Parents are equall!!");
                //System.out.println("parent1 "+parent1);
                //System.out.println("parent2 "+parent2);
                continue;            	
            }

            for (int j = 0; j < crossPos.size(); j++) {
                int cand1, cand2;
                int k = crossPos.get(j);
                
                if (rng.nextDouble() < 0.5D) {
                	mutateGeneCL(CL1, offspring1, k);
                	mutateGeneCL(CL2, offspring2, k);
                } 
                
                CL1 = updateCL(offspring1);
                CL2 = updateCL(offspring2);
            }

//            offspring1.calcFitness(ObjFunction);
//            offspring2.calcFitness(ObjFunction);

            if (offspringMinDist(offsprings, offspring1) == 0) {
            	// existe um individuo igual ao offspring1, mutate offspring1 and add
            	mutateGeneCL(CL1, offspring1, rng.nextInt(chromosomeSize));
            	offspring1.calcFitness(ObjFunction);
            	offsprings.add(offspring1);
            } else {
            	offspring1.calcFitness(ObjFunction);
            	offsprings.add(offspring1);
            }
            if (offspringMinDist(offsprings, offspring2) == 0) {
            	// existe um individuo igual ao offspring2, mutate offspring2 and add
            	mutateGeneCL(CL2, offspring2, rng.nextInt(chromosomeSize));
            	 offspring2.calcFitness(ObjFunction);
            	offsprings.add(offspring2);
            } else {
            	 offspring2.calcFitness(ObjFunction);
            	offsprings.add(offspring2);
            }

        }

        return offsprings;
    }

    
    protected void mutateGeneCL(ArrayList<Integer> CL, Chromosome<Integer> chromosome, Integer locus) {
        if (chromosome.get(locus) == 1) {
        	this.tripleElements[locus].qttUsed--;
        	this.decrementEsperanca();
            chromosome.set(locus, 0);
        } else if (CL.contains(locus)) {
        	this.tripleElements[locus].qttUsed++;
        	this.incrementEsperanca();
            chromosome.set(locus, 1);
        }
    }

    @Override
    public Boolean mutationCriteria() {
//    	calcular desvio padrão aqui
    	double desvioPadrao = 0d;
    	for (int i = 0; i < tripleElements.length; i++) {
    		int qtt = tripleElements[i].getQttUsed();
    		desvioPadrao += (qtt * qtt) - (2 * qtt * esperanca) + (esperanca * esperanca); 
		}
    	desvioPadrao = Math.sqrt(desvioPadrao/tripleElements.length);
    	
    	return rng.nextDouble() < mutationRate;
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
            for (int i = 0; i < chromosomeSize; i++) {
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
        GA_QBFPT ga = new GA_QBFPT(30, 1000, 100, 1.0 / 100.0, "instances/qbf020", AbstractGA.DEFAULT_CROSSOVER, AbstractGA.DEFAULT_MUTATION);
        Solution<Integer> bestSol = ga.solve();
        System.out.println("maxVal = " + bestSol);
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

    }

	@Override
	public Solution<Integer> createEmptySol() {
		Solution<Integer> sol = new Solution<Integer>();
        sol.cost = 0.0;
        return sol;
	}

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
		this.resetTripleElementsQttUsed();
		
	}

	@Override
	protected void mutateGene(Chromosome<Integer> chromosome, Integer locus) {
        if (chromosome.get(locus) == 1) {
        	this.tripleElements[locus].qttUsed--;
        	this.decrementEsperanca();
            chromosome.set(locus, 0);
        } else if (updateCL(chromosome).contains(locus)) {
        	this.tripleElements[locus].qttUsed++;
        	this.incrementEsperanca();
            chromosome.set(locus, 1);
        }    
		
	}

}
