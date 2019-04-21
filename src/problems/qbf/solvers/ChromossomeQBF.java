/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package problems.qbf.solvers;

import metaheuristics.ga.Chromosome;
import problems.Evaluator;
import solutions.Solution;

/**
 *
 * @author Cristina Bazzano
 * @author Jônatas Trabuco Belotti [jonatas.t.belotti@hotmail.com]
 * @author Matheus Diógenes Andrade
 */
public class ChromossomeQBF extends Chromosome<Integer> {

    @Override
    public void calcFitness(Evaluator<Integer> objEval) {
        Solution<Integer> solution = new Solution<Integer>();
        solution.cost = 0.0;

        for (int locus = 0; locus < size(); locus++) {
            if (get(locus) == 1) {
                solution.add(new Integer(locus));
            }
        }

        objEval.evaluate(solution);
        this.fitnessVal = solution.cost;
    }

}
