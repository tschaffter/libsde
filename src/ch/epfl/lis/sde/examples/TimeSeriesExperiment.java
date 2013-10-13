/*
Copyright (c) 2009-2012 Thomas Schaffter

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper(s) listed
on http://tschaffter.ch/projects/libsde.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

package ch.epfl.lis.sde.examples;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.net.URI;

import com.esotericsoftware.minlog.Log;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.epfl.lis.sde.SdeSettings;
import ch.epfl.lis.sde.solvers.SdeSolver;
import ch.epfl.lis.sde.solvers.SdeSolverFactory;

/**
 * Illustrates how to integrate a system of N stochastic differential equations (SDE) using libSDE.
 *
 * @version June 28, 2009
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class TimeSeriesExperiment {

	/** Time series data. */
	private DoubleMatrix2D timeSeries_;
	/** Time scale.*/
	private DoubleMatrix1D timeScale_;
	/** The duration of the experiment. */
	private double maxt_;
	/** Time interval between two time points (maxt/(numTimePoints-1)). */
	private double dt_;
	/** Number of time points. */
	private int numTimePoints_;
	
	/** Method used to integrate SDEs. */
	private SdeSolver solver_;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/**  */
	public static void main(String[] args) {
		
		@SuppressWarnings("unused")
		SdeSettings settings = SdeSettings.getInstance();
		
		Log.info("libSDE: Java library for integration of systems of stochastic differential equations (SDE).\n");
		
		try {
			int N = 50;
			Log.info("TimeSeriesExperiment", "Integrating N=" + N + " stochastic differential equations");
			TimeSeriesExperiment ts = new TimeSeriesExperiment();
			ts.run(N);
			
		} catch (OutOfMemoryError e) {
			String message = "There is not enough memory available to run this program.";
			message += "Quit one or more programs, and then try again.";
			message += "If enough amounts of RAM are installed on this computer, " +
				"try to run the Java Virtual Machine (JVM) with the command-line argument -Xmx1024m to use maximum 1024Mb of memory, " +
				"-Xmx2048m to use max 2048Mb, etc.";
			Log.error("TimeSeriesExperiment", message, e);
		} catch (Exception e) {
			Log.error("TimeSeriesExperiment", "Error while running the time series example.", e);
		}
		
		Log.info("");
		Log.info("Author: Thomas Schaffter (thomas.schaff...@gmail.com)");
		Log.info("Website: http://tschaffter.ch/projects/libsde");
		Log.info("");
		Log.info("Copyright (c) 2009-2012 Thomas Schaffter");
	}
	
	// ----------------------------------------------------------------------------
	
	/** Default constructor. */
	public TimeSeriesExperiment() {

		maxt_ = 0.;
		dt_ = 0.;
		numTimePoints_ = 0;
		timeSeries_ = null;
		timeScale_ = null;
		solver_ = null;
	}
	
	// ----------------------------------------------------------------------------
	
	/** Instantiates the different variables _before_ simulate the experiment. */
	public void init(int numTimePoints) {
		
		numTimePoints_ = numTimePoints;

		if (numTimePoints_ <= 1)
			throw new IllegalArgumentException("The number of time points must be greater than 1.");
		
		// allocation
		int n = solver_.getSystem().getDimension(); // dimension of the system
		timeSeries_ = new DenseDoubleMatrix2D(numTimePoints_, n);
		timeScale_ = new DenseDoubleMatrix1D(numTimePoints_);
		
		// initialize solver
		solver_.initialize();
	}
	
	// ----------------------------------------------------------------------------

	/**
	 * Run the experiment.
	 * @param N Number of stochastic differential equations to solve. 
	 */
	public void run(int N) {

		// set the parameters first
		SdeSettings settings = SdeSettings.getInstance();
		settings.setSeed(-1); // random seed
		settings.setMaxt(1);
		settings.setDt(0.001);
		
		// set SDEs
		ExampleSde system = new ExampleSde(N);
		system.setName("mySDE");
		system.setSigma(0.2);
		
		// set solver
		solver_ = SdeSolverFactory.createSolver(SdeSolverFactory.SRK_ITO);
		if (solver_ == null)
			throw new RuntimeException("Unable to instantiate the solver");
		
		DoubleMatrix1D X0 = new DenseDoubleMatrix1D(system.getDimension()); // initial condition
		solver_.setX(X0.assign(1)); // set initial condition X0=[1,1,1,...]^T
		
		// associate the system to the solver
		solver_.setSystem(system);
		Log.info("TimeSeriesExperiment", solver_.getDescription() + "\n");
		

		// set the time series experiment
		init(101); // number of time points wished
		
		integrate();		
		
		// save time series to file
		URI target = new File(system.getName() + "_timeseries.tsv").toURI();
		Log.info("TimeSeriesExperiment", "Writing time series to " + target.getPath());
		printAll(target);
	}
		
	// ----------------------------------------------------------------------------
	
	/** Numerical integration of the system of N SDEs.*/
	public void integrate() {
		
		SdeSettings settings = SdeSettings.getInstance();
		
		int n = solver_.getSystem().getDimension();
		double t = 0;
		maxt_ = settings.getMaxt();
		
		if (maxt_ <= 0)
			throw new IllegalArgumentException("Duration (maxt) must be greater than 0.");
		
		dt_ = maxt_/(double)(numTimePoints_-1);
		
		if (dt_ <= 0 || dt_ > maxt_)
			throw new IllegalArgumentException("Interval between two measured time points must be greater than 0 and and smaller than maxt.");
		
		double frac = dt_/settings.getDt();
		if (frac - (int)frac != 0) // check whether the decimal part of frac is 0 or not
			throw new IllegalArgumentException("Interval between two measured time points (maxt/(numTimePoints-1)) must be a multiple of the integration step-size.");
		
		solver_.setH(dt_);
		
		DoubleMatrix1D X = solver_.getX();
		if (X == null)
			throw new RuntimeException("TimeSeriesExperiment.integrate(): X0 null");
		
		// Set first line of the time series dataset (at t=0)
		for (int i=0; i<n; i++)
			timeSeries_.set(0, i, X.get(i));
		
		int pt = 1;
		
		do {
			double t1 = t;
			try {

				// this steps the time by TimeSeriesExperiment.dt_, the solver integrates with a smaller, fixed step size
				// defined in SdeSettings by dt_*multiplier_ (SdeSettings.dt_ != TimeSeriesExperiment.dt_)
				// WARNING: TimeSeriesExperiment.dt_ must be a multiple of SdeSettings.dt_
				t += solver_.step();

			} catch (Exception e) {
				throw new RuntimeException("TimeSeriesExperiment.integrate(): Exception at t = " + t + ": " + e.getMessage());
			}
			
			if (t != t1 + dt_)
				throw new RuntimeException("TimeSeriesExperiment.integrate(): Solver failed to step time by dt, expected t = " + (t1+dt_) + ", obtained t = " + t);

			// save the result
			X = solver_.getX();
			for (int i = 0; i < n; i++)
				timeSeries_.set(pt, i,X.get(i));
			timeScale_.set(pt, t);
			
			pt++;
			
		} while (t < maxt_);

		assert t == maxt_ : "t=" + t + " maxt=" + maxt_;
		assert pt == numTimePoints_;
	}
	
	// ----------------------------------------------------------------------------

	/** Wrapper function to print trajectories to files. */
	public void printAll(URI target) {
		
		try {
			printTrajectories(target, timeSeries_, timeScale_);
		} catch (Exception e) {
			Log.error("TimeSeriesExperiment", "Failed to write time series to file.", e);
		}
	}
	
	// ----------------------------------------------------------------------------

	/** Writes time series data with the time scale to file. */
	public void printTrajectories(URI target, DoubleMatrix2D timeSeries, DoubleMatrix1D timeScale) throws IOException, Exception {

			FileWriter fw = new FileWriter(target.getPath(), false);
			int R = timeSeries.rows();
			int C = timeSeries.columns();
			
			for (int i = 0; i < R; i++) {
				// first column of the file is time scale
				fw.write(Double.toString(timeScale.get(i)));
				
				for (int j = 0; j < C; j++)
					fw.write("\t" + timeSeries.get(i, j));
				
				fw.write("\n");
			}
			fw.close();
	}
	
	// ============================================================================
	// GETTERS AND SETTERS

	public void setNumTimePoints(int numTimePoints) { numTimePoints_ = numTimePoints; }
	public int getNumTimePoints() { return numTimePoints_; }
	
	public DoubleMatrix2D getTimeSeries() { return timeSeries_; }
	
	public void setSolver(SdeSolver solver) { solver_ = solver; }
	public SdeSolver getSolver() { return solver_; }
	
	public void setMaxt(double maxt) { maxt_ = maxt; }
	public double getMaxt() { return maxt_; }
}
