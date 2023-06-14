import java.util.*;

public class ExampleExperiment {

	/**
	 * The maximal budget for evaluations done by an optimization algorithm equals
	 * dimension * BUDGET_MULTIPLIER.
	 * Increase the budget multiplier value gradually to see how it affects the runtime.
	 */
	public static final int BUDGET_MULTIPLIER = 4;

	/**
	 * The maximal number of independent restarts allowed for an algorithm that restarts itself.
	 */
	public static final int INDEPENDENT_RESTARTS = 10000;

	/**
	 * The random seed. Change if needed.
	 */
	public static final long RANDOM_SEED = 123;//0xdeadbeef;

	/**
	 * The problem to be optimized (needed in order to simplify the interface between the optimization
	 * algorithm and the COCO platform).
	 */
	public static Problem PROBLEM;
	/**
	 * Interface for function evaluation.
	 */
	public interface Function {
		double[] evaluate(double[] x);
		double[] evaluateConstraint(double[] x);
	}

	/**
	 * Evaluate the static PROBLEM.
	 */
	public static final Function evaluateFunction = new Function() {
		public double[] evaluate(double[] x) {
			return PROBLEM.evaluateFunction(x);
		}
		public double[] evaluateConstraint(double[] x) {
			return PROBLEM.evaluateConstraint(x);
		}
	};
	/**
	 * The main method initializes the random number generator and calls the example experiment on the
	 * bi-objective suite.
	 */
	public static void main(String[] args) {

		Random randomGenerator = new Random(RANDOM_SEED);

		/* Change the log level to "warning" to get less output */
		CocoJNI.cocoSetLogLevel("info");

		System.out.println("Running the example experiment... (might take time, be patient)");
		System.out.flush();

		/* Start the actual experiments on a test suite and use a matching logger, for
		 * example one of the following:
		 *
		 *   bbob                 24 unconstrained noiseless single-objective functions
		 *   bbob-biobj           55 unconstrained noiseless bi-objective functions
		 *   [bbob-biobj-ext       92 unconstrained noiseless bi-objective functions]
		 *   bbob-largescale      24 unconstrained noiseless single-objective functions in large dimension
		 *   [bbob-constrained*   48 constrained noiseless single-objective functions]
		 *   bbob-mixint          24 unconstrained noiseless single-objective functions with mixed-integer variables
		 *   bbob-biobj-mixint    92 unconstrained noiseless bi-objective functions with mixed-integer variables
		 *
		 * Suites with a star are partly implemented but not yet fully supported.
		 *
		 * Adapt to your need. Note that the experiment is run according
		 * to the settings, defined in exampleExperiment(...) below.
		 */
		dGExperiment("bbob", "bbob", randomGenerator);

		System.out.println("Done!");
		System.out.flush();

		return;
	}
	/**
	 * A simple example of benchmarking random search on a given suite with default instances
	 * that can serve also as a timing experiment.
	 *
	 * @param suiteName Name of the suite (e.g. "bbob", "bbob-biobj", or "bbob-constrained").
	 * @param observerName Name of the observer matching with the chosen suite (e.g. "bbob-biobj"
	 * when using the "bbob-biobj-ext" suite).
	 * @param randomGenerator The random number generator.
	 */
	public static void dGExperiment(String suiteName, String observerName, Random randomGenerator) {
		try {

			/* Set some options for the observer. See documentation for other options. */
			final String observerOptions =
					"result_folder: RS_on_" + suiteName + " "
							+ "algorithm_name: RS "
							+ "algorithm_info: \"A simple random search algorithm\"";

			/* Initialize the suite and observer.
			 * For more details on how to change the default options, see
			 * http://numbbo.github.io/coco-doc/C/#suite-parameters and
			 * http://numbbo.github.io/coco-doc/C/#observer-parameters. */
			Suite suite = new Suite(suiteName, "", "");
			Observer observer = new Observer(observerName, observerOptions);
			Benchmark benchmark = new Benchmark(suite, observer);

			/* Initialize timing */
			Timing timing = new Timing();

			/* Iterate over all problems in the suite */
			while ((PROBLEM = benchmark.getNextProblem()) != null) {

				int dimension = PROBLEM.getDimension();

				/* Run the algorithm at least once */
				for (int run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

					long evaluationsDone = PROBLEM.getEvaluations() + PROBLEM.getEvaluationsConstraints();
					long evaluationsRemaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluationsDone;

					/* Break the loop if the target was hit or there are no more remaining evaluations */
					if (PROBLEM.isFinalTargetHit() || (evaluationsRemaining <= 0))
						break;

					/* Call the optimization algorithm for the remaining number of evaluations */
					int populationSize = 1000;
					double crossoverRate = 0.5;
					double factor = randomGenerator.nextDouble() * 2.0;
					algorithm(evaluateFunction,
							dimension,
							PROBLEM.getNumberOfObjectives(),
							PROBLEM.getNumberOfConstraints(),
							PROBLEM.getSmallestValuesOfInterest(),
							PROBLEM.getLargestValuesOfInterest(),
							PROBLEM.getNumberOfIntegerVariabls(),
							evaluationsRemaining,
							randomGenerator,
							populationSize,
							factor,
							crossoverRate);

					/* Break the loop if the algorithm performed no evaluations or an unexpected thing happened */
					if (PROBLEM.getEvaluations() == evaluationsDone) {
						System.out.println("WARNING: Budget has not been exhausted (" + evaluationsDone + "/"
								+ dimension * BUDGET_MULTIPLIER + " evaluations done)!\n");
						break;
					} else if (PROBLEM.getEvaluations() < evaluationsDone)
						System.out.println("ERROR: Something unexpected happened - function evaluations were decreased!");
				}

				/* Keep track of time */
				timing.timeProblem(PROBLEM);
			}

			/* Output the timing data */
			timing.output();

			benchmark.finalizeBenchmark();

		} catch (Exception e) {
			System.err.println(e.toString());
		}
	}

	/**
	 * Evolution algorithm.
	 */
	public static void algorithm(Function f,
								 int dimension,
								 int numberOfObjectives,
								 int numberOfConstraints,
								 double[] lowerBounds,
								 double[] upperBounds,
								 int numberOfIntegerVariables,
								 long maxBudget,
								 Random randomGenerator,
								 int populationSize,
								 double factor,
								 double crossoverRate) {

		// Set dLow and dHigh, values according to literature
		double dLow = 0.25;
		double dHigh = 0.25;
		boolean exploit = true; // false means explore

		// Initialize population
		double[][] population = new double[populationSize][dimension];
		for (int i = 0; i < populationSize; i++) {
			population[i] = generateRandomSolution(lowerBounds, upperBounds, numberOfIntegerVariables, randomGenerator);
		}

		// Initialize history
		List<double[]> history = new ArrayList<>();
		history.addAll(Arrays.asList(population));

		// Main loop
		int t = 0;
		while (t < maxBudget) {

			// Calculate diversity and set mode
			double diversity = diversity(population, populationSize, dimension, lowerBounds, upperBounds);
			if (diversity <= dLow)
				exploit = false;
			else if (diversity > dHigh)
				exploit = true;

			// For each individual in the population
			for (int i = 0; i < populationSize; i++) {

				// Select a random individual as the base vector
				//double[] baseVector = population[randomGenerator.nextInt(populationSize)];
				//double[] targetVector = population[i];

				if (exploit)
				{
					// Selection (above if-else) and recombination with other vector from population
					// Select a random individual as the base vector
					double[] baseVector = population[randomGenerator.nextInt(populationSize)];

					// Sample two random individuals as the difference vectors
					double[] differenceVector1 = population[randomGenerator.nextInt(populationSize)];
					//double[] differenceVector2 = population[randomGenerator.nextInt(populationSize)];

					// Generate a mutant vector by adding the scaled difference vectors to the base vector
					//double[] mutantVector = new double[dimension];
					//for (int j = 0; j < dimension; j++) {
					//	mutantVector[j] = baseVector[j] + factor * (differenceVector1[j] - differenceVector2[j]);
					//}

					// Perform crossover between the mutant vector and the target vector
					double[] targetVector = population[i];
					double[] trialVector = BinaryCrossover(differenceVector1, baseVector, crossoverRate, randomGenerator);

					// Update the population and history
					history.add(trialVector);
					if (numberOfConstraints>0)
					{
						f.evaluateConstraint(targetVector);
						f.evaluateConstraint(trialVector);

					}
					population[i] = Tournament(targetVector, trialVector, f); // brak tournament powoduje "budget has not been exhausted"
				}
				else
				{
					// Mutation
					// Select a random individual as the base vector
					double[] baseVector = population[randomGenerator.nextInt(populationSize)];

					// Sample two random individuals as the difference vectors
					double[] differenceVector1 = population[randomGenerator.nextInt(populationSize)];
					double[] differenceVector2 = population[randomGenerator.nextInt(populationSize)];

					// Generate a mutant vector by adding the scaled difference vectors to the base vector
					double[] mutantVector = new double[dimension];
					for (int j = 0; j < dimension; j++) {
						mutantVector[j] = baseVector[j] + factor * (differenceVector1[j] - differenceVector2[j]);
					}

					// Perform crossover between the mutant vector and the target vector
					double[] targetVector = population[i];
					//double[] trialVector = BinaryCrossover(targetVector, mutantVector, crossoverRate, randomGenerator);

					// Update the population and history
					//history.add(trialVector);
					history.add(mutantVector);
					if(numberOfConstraints>0)
					{
						f.evaluateConstraint(mutantVector);
					}
					population[i] = mutantVector; //Tournament(targetVector, mutantVector, f); // lub bez Tournament?
				}
			}

			t++;
		}
	}

	public static double[] BinaryCrossover(double[] x, double[] y, double cr, Random randomGenerator)
	{
		double[] z = new double[x.length];
		double a = 0;

		for(int i = 0; i < x.length; i++)
		{
			a = randomGenerator.nextDouble();
			if(a < cr)
				z[i] = y[i];
			else
				z[i] = x[i];
		}

		return z;
	}

	public static double[] Tournament(double[] x, double[] y, Function fitness)
	{
		// The objective vector that is the result of the evaluation (in single-objective problems only the first vector item is being set).
		if(fitness.evaluate(x)[0] < fitness.evaluate(y)[0])
			return x;
		else
			return y;
	}

	public static double[] generateRandomSolution(double[] lowerBounds, double[] upperBounds, int numberOfIntegerVariables, Random randomGenerator) {

		int dimension = lowerBounds.length;
		double[] solution = new double[dimension];

		for (int i = 0; i < dimension; i++) {
			if (i < numberOfIntegerVariables) {
				// Generate a random integer value for integer variables
				solution[i] = randomGenerator.nextInt((int) (upperBounds[i] - lowerBounds[i] + 1)) + lowerBounds[i];
			} else {
				// Generate a random real value for real variables
				solution[i] = lowerBounds[i] + randomGenerator.nextDouble() * (upperBounds[i] - lowerBounds[i]);
			}
		}

		return solution;
	}

	public static double diversity(double[][] population, int populationSize, int dimensionality,
								   double[] lowerBounds, double[] upperBounds)
	{
		double diversity = 0;
		double sum = 0;
		double[] averagePoint = averagePoint(population, populationSize, dimensionality);
		for (int i = 0; i < populationSize; i++)
		{
			double sumOfValues = 0;
			for (int j = 0; j < dimensionality; j++)
			{
				sumOfValues += Math.pow((population[i][j] - averagePoint[j]), 2);
			}
			sum += Math.sqrt(sumOfValues);
		}
		diversity = sum / (populationSize * diagonalOfSearchSpace(lowerBounds, upperBounds));

		return diversity;
	}

	public static double[] averagePoint(double[][] population, int populationSize, int dimensionality)
	{
		double[] averagePoint = new double[dimensionality];

		for (int i = 0; i < dimensionality; i++)
		{
			double averageValue = 0;
			for (int j = 0; j < populationSize; j++)
			{
				averageValue += population[j][i];
			}
			averageValue /= populationSize;
			averagePoint[i] = averageValue;
		}

		return averagePoint;
	}

	public static double diagonalOfSearchSpace(double[] lowerBounds, double[] upperBounds)
	{
		double L = 0;
		for (int i = 0; i < lowerBounds.length; i++)
		{
			L += Math.pow((lowerBounds[i] - upperBounds[i]),2);
		}
		return Math.sqrt(L);
	}
}