/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package problems.qbfpt.solvers;

import java.io.IOException;
import problems.qbf.solvers.GA_QBF;

/**
 *
 * @author
 */
public class GA_QBFPT extends GA_QBF {

    public GA_QBFPT(Integer tempoExecucao, Integer geracoesConvengencia, Integer popSize, Double mutationRate, String filename) throws IOException {
        super(tempoExecucao, geracoesConvengencia, popSize, mutationRate, filename);
    }

}
