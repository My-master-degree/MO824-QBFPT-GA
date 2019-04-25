package main;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import problems.qbf.QBF;
import problems.qbfpt.solvers.GA_QBFPT;
import solutions.Solution;

/**
 * Class that executes the GA for MAXQBFPT problem
 *
 *
 * @author Cristina Bazzano
 * @author Jônatas Trabuco Belotti [jonatas.t.belotti@hotmail.com]
 * @author Matheus Diógenes Andrade
 */
public class Main {

    public static int TIME_LIMIT = 30; // In minutes
    public static int GENERATIONS_LIMIT = -1; // Generations without improvement in incumbent (negative values for not using iterations limit)
    public static String outputCsv;

    public static final String[] FILES_LIST = new String[]{
        "instances/qbf020",
        "instances/qbf040",
        "instances/qbf060",
        "instances/qbf080",
        "instances/qbf100",
        "instances/qbf200",
        "instances/qbf400"
    };

    //Calls execution method with 5 different configurations
    public static void main(String[] args) throws IOException {

        outputCsv = "fileName,tamanhoPop,taxaMut,tipoCrossover,tipoMutacao,manutencaoDiversidade, solucao\n";

        executeGA(1, 1D/100D, GA_QBFPT.XOR_UNIFORM_CROSSOVER, GA_QBFPT.DEFAULT_MUTATION, false); // CONF P -> Padrão
        executeGA(2, 1D/100D, GA_QBFPT.XOR_UNIFORM_CROSSOVER, GA_QBFPT.DEFAULT_MUTATION, false); // CONF A -> Padrão + tamanho população alternativo
        executeGA(1, 5D/100D, GA_QBFPT.XOR_UNIFORM_CROSSOVER, GA_QBFPT.DEFAULT_MUTATION, false); // CONF B -> Padrão + taxa mutação alternativa
        executeGA(1, 1D/100D, GA_QBFPT.XOR_CROSSOVER, GA_QBFPT.DEFAULT_MUTATION, false); // CONF C -> Padrão + crossover XOR
        executeGA(1, 1D/100D, GA_QBFPT.DEFAULT_CROSSOVER, GA_QBFPT.DEFAULT_MUTATION, false); // CONF D -> Padrão + crossover de 2 pontos
        executeGA(1, 1D/100D, GA_QBFPT.XOR_UNIFORM_CROSSOVER, GA_QBFPT.DYNAMIC_MUTATION, false); // CONF E -> Padrão + Mutação adaptativa
        executeGA(1, 1D/100D, GA_QBFPT.XOR_UNIFORM_CROSSOVER, GA_QBFPT.DEFAULT_MUTATION, true); // CONF F -> Padrão + Manutenção da diversidade

        saveOutput("output.csv", outputCsv);
    }

    private static void executeGA(int tamPop, double taxaMuta, int tipoCrossover, int tipoMutacao, boolean manutencaoDiversidade) throws IOException {

        long beginTotalTime = System.currentTimeMillis();

        // Iterating over files
        for (String arquivo : FILES_LIST) {

            //Print configurations of the execution
            System.out.println("Executing GA for file: " + arquivo);
            System.out.println("Configuration:");

            printStopCriterion();

            // Executing GA
            System.out.println("Execution:");

            long beginInstanceTime = System.currentTimeMillis();
            
            GA_QBFPT ga = new GA_QBFPT(TIME_LIMIT, GENERATIONS_LIMIT, new QBF(arquivo).getDomainSize() * tamPop, taxaMuta, arquivo, tipoCrossover, tipoMutacao, manutencaoDiversidade);
            Solution<Integer> bestSolution = ga.solve();
            System.out.println(" maxVal = " + bestSolution);

            long endInstanceTime = System.currentTimeMillis();
            long totalInstanceTime = endInstanceTime - beginInstanceTime;
            System.out.println("Time = " + (double) totalInstanceTime / (double) 1000 + " seg");
            System.out.println("\n");

            // "fileName,tamanhoPop,taxaMut,tipoCrossover,tipoMutacao,manutencaoDiversidade,solucao
            outputCsv += arquivo + "," + tamPop + "," + taxaMuta + "," + tipoCrossover + "," + tipoMutacao + "," + manutencaoDiversidade + "," + bestSolution.cost + "\n";

        }

        // Calculating time of all executions
        long totalTime = System.currentTimeMillis() - beginTotalTime;

        System.out.println("Tempo execução todos arquivos: " + (totalTime / 1000D) + "seg \n"
                + "----------------------------------------------------- \n \n");
    }

    private static void printStopCriterion() {
        String resp = " Stop Criterion = ";

        resp += TIME_LIMIT + " minutes";

        if (GENERATIONS_LIMIT > 0) {
            resp += " ou " + GENERATIONS_LIMIT + " generations without new incumbent";
        }

        System.out.println(resp);
    }

    public static void saveOutput(String fileName, String content) {
        File dir;
        PrintWriter out;

        dir = new File("output");

        if (!dir.exists()) {
            dir.mkdirs();
        }

        try {
            out = new PrintWriter(new File(dir, fileName));
            out.print(content);
            out.close();
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }
    }
}
