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

package ch.epfl.lis.sde.solvers;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.random.Normal;
import ch.epfl.lis.sde.Sde;
import ch.epfl.lis.sde.SdeSettings;

/** 
 * This class is extended by the SDE solvers implemented in libSDE.
 * 
 * @version June 28, 2009
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public abstract class SdeSolver {

	/** SDEs system. */
	protected Sde system_;
	
	/** Wiener increments dW. */
	protected DoubleMatrix2D W_;
	/** dZ increments. */
	protected DoubleMatrix2D Z_;
	/** Approximation of the solution of the SDEs system. */
	protected DoubleMatrix1D X_;

	/** Current time. */
	protected double currentTime_;
	/** Number of time points for the Wiener process. */ // TODO: In one interval ?
	protected int numTimePointsWiener_;
	/** Drift coefficients. */
	protected DoubleMatrix1D F_;
	/** Diffusion coefficients. */
	protected DoubleMatrix2D G_;
	/** Internal integration step size. */
	protected double h_;
	/** External integration step  size. */
	protected double H_;
	
	/** Absolute _or_ relative precision _per variable_ need to be satisfied for convergence. */
	private double absolutePrecision_ = 0.000001;		// Not yet implemented
	/** See absolutePrecision_, in addition, this is also the tolerance used for integration. */ 
	private double relativePrecision_ = 0.0001;			// Not yet implemented
	
	/** Current number of evaluations of the system. */
	private int numEvaluations_;
	
	/** Is true if the system has converged. */
	private boolean converged_;							// Not yet implemented	

	// ============================================================================
	// PUBLIC METHODS
	
	/** Proceed on iteration of the numerical integration (step size = h). */
	public abstract void advance(final double t, final double h, final DoubleMatrix1D dW, 
			final DoubleMatrix1D dZ, final DoubleMatrix1D Xin, DoubleMatrix1D Xout) throws Exception;
	
	// ----------------------------------------------------------------------------
	
	/** Default constructor. */
	public SdeSolver() {
		
		system_ = null;
		reset();
	}
	
	// ----------------------------------------------------------------------------
	
	/** Reset the solver. */
	public void reset() {
		
		W_ = null;
		Z_ = null;
		X_ = null;
		currentTime_ = 0;
		numTimePointsWiener_ = 0;
		F_ = null;
		G_ = null;
		h_ = 0.;
		H_ = 0.;
		numEvaluations_ = 0;
	}
	
	// ----------------------------------------------------------------------------
	
	/** At each step(), W must be regenerated. */
	public void generateW() {
		
		Normal normal = SdeSettings.getInstance().getNormalDistribution();
		double sqrt_dt = Math.sqrt(SdeSettings.getInstance().getDt());
		
		for (int i = 0; i < numTimePointsWiener_; i++) {
			for (int j = 0; j < system_.getDimension(); j++) {
				double N1 = normal.nextDouble();
				W_.set(i, j, N1*sqrt_dt);
			}
		}
	}

	// ----------------------------------------------------------------------------
	
	/** At each step(), regenerate W and Z (Z are useful for Runge-Kutta of strong order >= 1.5). */
	public void generateWZ() {
		
		SdeSettings settings = SdeSettings.getInstance();
		double sqrt_dt = Math.sqrt(settings.getDt());
		double sqrt_3 = Math.sqrt(3.);
		double pow_dt_3_2 = Math.pow(settings.getDt(), 3.0/2.0);
		double N1, N2;

		for (int i = 0; i < numTimePointsWiener_; i++) {
			for (int j = 0; j < system_.getDimension(); j++) {
				N1 = settings.getNormalDistribution().nextDouble();
				N2 = settings.getNormalDistribution().nextDouble();
				W_.set(i, j, N1*sqrt_dt);
				Z_.set(i, j, 0.5*(N1+1/sqrt_3*N2)*pow_dt_3_2);
			}
		}
	}
	
	// ----------------------------------------------------------------------------
	
	/**
	 * <p>Must be call after the SDE parameterized and associated to this solver.<br>
	 * The correct seed and all the SdeSettings parameters must also have been set.</p>
	 */
	public void initialize() {
		
		SdeSettings settings = SdeSettings.getInstance();
		
		if (system_ == null)
			throw new IllegalArgumentException("No system of SDEs has been set!");
		
		if (settings.getMultiplier() == 0 || settings.getDt() == 0)
			throw new RuntimeException("Nothing to do!");
	
		// Number of Wiener increments in one integration step-size h_
		numTimePointsWiener_ = settings.getMultiplier();
		h_ = settings.getDt()*(double)settings.getMultiplier();
		currentTime_ = 0; // integration start from 0
		
		int n = system_.getDimension();
		W_ = new DenseDoubleMatrix2D(numTimePointsWiener_, n);
		Z_ = new DenseDoubleMatrix2D(numTimePointsWiener_, n);
		
		settings.initializeRNG(); // use the defined seed to set up the RNG
		
		// drift vector and diffusion matrix
		F_ = new DenseDoubleMatrix1D(n);
		G_ = new DenseDoubleMatrix2D(n, n);
	}
	
	// ----------------------------------------------------------------------------
	
	/** Step the integration from the current time t1 to t1+H_, return H_. */
	public double step() throws Exception {
		
		double t1 = currentTime_; // get the current time
		double t1_bkp = t1; // to compute the effective step size in the end
		double t2 = t1 + H_; // get the time at which step() must return

		int n = system_.getDimension();
		DoubleMatrix1D Xout = new DenseDoubleMatrix1D(n);
		DoubleMatrix1D dW = new DenseDoubleMatrix1D(n);
		DoubleMatrix1D dZ = new DenseDoubleMatrix1D(n);
		
		while (t1 < t2) {
			
			//System.out.println("t1: " + t1);
			generateWienerIncrements(n, dW, dZ); // generate new independent Wiener process samples
			
			// compute the drift and the diffusion
			system_.getDriftAndDiffusion(t1, X_, F_, G_);
			
			// compute the next approximation
			advance(t1, h_, dW, dZ, X_, Xout);
			
			// if required, perform a check on X
			checkX(Xout);
			
			// checks if the system converged
			// NOT YET IMPLEMENTED but the function checkConvergence(Xout) can be ower-written by the user
			checkConvergence(Xout);
			
			// save the current solution
			X_.assign(Xout);
			
			t1 += h_;
			numEvaluations_++;
		}
		
		currentTime_ = t2;
		return currentTime_ - t1_bkp;
	}
	
	// ----------------------------------------------------------------------------
	
	/**
	 * <p>Generates Wiener increments.<br>
	 * <b>Important:</b> Override this function if dW is required (e.g. SRK methods)</p>
	 */
	public void generateWienerIncrements(int n, DoubleMatrix1D dW, DoubleMatrix1D dZ) {
		
		generateW();
		
		for (int j = 0; j < n; j++) {
			dW.set(j, 0);
			for (int k = 0; k < W_.rows(); k++)
				dW.set(j, dW.get(j)+W_.get(k, j));
		}
	}
	
	// ----------------------------------------------------------------------------
	
	/** This function can be overwritten if a check should be performed on X at every step(). */
	public void checkX(DoubleMatrix1D X) {}
	
	// ----------------------------------------------------------------------------
	
	/** This function checks the convergence of the system. If yes, set converged_ = true. */
	public void checkConvergence(DoubleMatrix1D X) {}
	
	// ----------------------------------------------------------------------------
	
	/**
	 * <p>SDE integration. Returns the solution X at t=maxt.</p>
	 * <br>
	 * <b>Important:</b> If you want to save the state of X at different time points between
	 * t=0 and t=maxt, you must redefine this function.
	 * 
	 * @see ch.epfl.lis.sde.examples.TimeSeriesExperiment
	 */
	public DoubleMatrix1D integrate() throws Exception {
		
		SdeSettings settings = SdeSettings.getInstance();

		double maxt = settings.getMaxt();
		H_ = maxt;
		
		initialize();
		
		//int pt = 1;
		
		do {
			try {
				step();
				//pt++;
				
			} catch (Exception e) {
				throw new Exception("SdeSolver():integrate(): Exception, t = " + currentTime_ + ": " + e.getMessage());
			}
		} while (currentTime_ < maxt);
		
		return X_;
	}
	
	// ----------------------------------------------------------------------------
	
//	/** 
//	 * Implements the method gsl_multiroot_test_delta() of GSL:
//	 * This function tests for the convergence of the sequence by comparing the last step dx
//	 * with the absolute error epsabs and relative error epsrel to the current position x.
//	 * The test returns true if the following condition is achieved:
//     * 		|dx_i| < epsabs + epsrel |x_i|
//     * for each component of x and returns false otherwise.
//	 */
//	public boolean converged(final DoubleMatrix1D X, final DoubleMatrix1D previousX) {
//		
//		int n = system_.getDimension();
//		double dx = 0.;
//		
//		for (int i=0; i<n; i++) {
//			
//			dx = Math.abs(X.get(i) - previousX.get(i));
//			
//			if (dx > absolutePrecision_ + relativePrecision_*Math.abs(X.get(i))) {				
//				return false;
//			}
//		}
//		return true;
//	}

	
    // =======================================================================================
    // GETTERS AND SETTERS
	
	public void setSystem(Sde sde) { system_ = sde; }
	public Sde getSystem() { return system_; }
	
	public void setX(DoubleMatrix1D X) { X_ = X; }
	public DoubleMatrix1D getX() { return X_; }
	
	public void setH(double H) { H_ = H;}
	public double getH() { return H_; }
	
	public void setAbsolutePrecision(double value) { absolutePrecision_ = value; }
	public double getAbsolutePrecision() { return absolutePrecision_; }
	
	public void setRelativePrecision(double value) { relativePrecision_ = value; }
	public double getRelativePrecision() { return relativePrecision_; }
	
	public boolean converged() { return converged_; }
	
	public String getDescription() { return "No description for this solver"; }
	
	public int getNumEvaluations() { return numEvaluations_; }
}
