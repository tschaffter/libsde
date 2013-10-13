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

/**
 * Implements the explicit Milstein method (strong order of convergence 1)to use with Stratonovich scheme.
 * 
 * <p>The iterative Milstein scheme is given by the formula described in [1].</p>
 * 
 * <p>Xn+1 = Xn + Fn*h + Gn*dWn + 0.5*Gn*Gn'*(dWn^2 - h) (Ito)<br>
 * Xn+1 = Xn + Fn*h + Gn*dWn + 0.5*Gn*Gn'*(dWn^2)     (Stratonovich)	<-- Implemented here</p>
 * 
 * <p>where h is the integration step size.</p>
 * 
 * <p>This solver implements a first-order stochastic Runge-Kutta algorithm (FSRK) that does
 * not need the user to specify explicitly the derivative of the diffusion term G' [1].</p>
 * 
 * <p>[1] P.E. Kloeden, E. Platen, and H. Schurz, Numerical solution of SDE through
 *     computer experiments, Springer, 1994. (pp 150-153)</p>
 * <p>[2] http://en.wikipedia.org/wiki/Runge-Kutta_method_(SDE)</p>
 *
 * @version June 28, 2009
 *
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class MilsteinStratonovich extends SdeSolver {

	/** Support vector. */
	private DoubleMatrix1D supportX_;
	/** Drift evaluated at the support vector. */
	private DoubleMatrix1D supportF_;
	/** Diffusion evaluated at the support vector. */
	private DoubleMatrix2D supportG_;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor. */
	public MilsteinStratonovich() {
		
		super();
	}
	
	// ----------------------------------------------------------------------------
	
	/** Reset the solver. */
	public void reset() {
		
		super.reset();
		supportX_ = null;
		supportF_ = null;
		supportG_ = null;
	}
	
	// ----------------------------------------------------------------------------
	
	/** Initialization of the solver. */
	public void initialize() throws IllegalArgumentException {
		
		super.initialize();
		int n = system_.getDimension();
		supportX_ = new DenseDoubleMatrix1D(n);
		supportF_ = new DenseDoubleMatrix1D(n);
		supportG_ = new DenseDoubleMatrix2D(n, n);
	}
	
	// ----------------------------------------------------------------------------
	
	/** Proceed on iteration of the numerical integration (step size = h). */
	public void advance(final double t, final double h, final DoubleMatrix1D dW,
			final DoubleMatrix1D dZ, final DoubleMatrix1D Xin, DoubleMatrix1D Xout) throws Exception {
		
		int n = system_.getDimension();
		double xout = 0.;
		double sqrt_h = Math.sqrt(h);
				
		// vector supporting value
		for (int i = 0; i < n; i++)
			supportX_.set(i, Xin.get(i) + F_.get(i)*h + G_.get(i, i)*sqrt_h);
		
		system_.getDriftAndDiffusion(t, supportX_, supportF_, supportG_);
		
		// for each dimension, we iterate the solution by h
		for (int i = 0; i < n; i++) {	
			xout = Xin.get(i) + F_.get(i)*h + G_.get(i, i)*dW.get(i) + (1.0/(2.0*sqrt_h))
				*(supportG_.get(i, i) - G_.get(i, i))*Math.pow(dW.get(i), 2.0);
			Xout.set(i, xout);
		}
	}
	
	// ----------------------------------------------------------------------------
	
	/** Returns the description of this solver. */
	public String getDescription() {
		
		String desc = "Explicit Milstein method for Stratonovich SDEs (Kloeden et al., 1994).\n";
		desc += "Strong order of convergence: 1\n";
		desc += "Weak order of convergence: 1\n";
		desc += "Only diagonal noise is handled. Avoids specification of derivatives (Runge-Kutta approach).";
		
		return desc;
	}
}
