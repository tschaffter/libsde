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

import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/** 
 * Implements the explicit Runge-Kutta method (strong order of convergence 1.5) to use with Ito scheme.
 * 
 * <p>Please refer to [1] p. 390</p>
 * 
 * <p>This solver implements a stochastic Runge-Kutta algorithm (SRK) that does
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
public class SRK15 extends SdeSolver {

	/** Require 4 supporting vectors. */
	private ArrayList<Support> s1_, s2_, s3_, s4_;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor. */
	public SRK15() {
		
		super();
	}
	
	// ----------------------------------------------------------------------------
	
	@Override
	public void generateWienerIncrements(int n, DoubleMatrix1D dW, DoubleMatrix1D dZ) {
		
		generateWZ();
		
		for (int j = 0; j < n; j++) {
			dW.set(j, 0);
			dZ.set(j, 0);
			for (int k = 0; k < W_.rows(); k++) { // dZ has the same size than dW
				dW.set(j, dW.get(j)+W_.get(k, j));
				dZ.set(j, dZ.get(j)+Z_.get(k, j));
			}
		}
	}
	
	// ----------------------------------------------------------------------------
	
	/** Reset the solver. */
	public void reset() {
		
		super.reset();
		s1_ = null;
		s2_ = null;
		s3_ = null;
		s4_ = null;
	}
	
	// ----------------------------------------------------------------------------
	
	/** Initialization of the solver. */
	public void initialize() throws IllegalArgumentException {
		
		super.initialize();
		int n = system_.getDimension();
		
		s1_ = new ArrayList<Support>(n);
		s2_ = new ArrayList<Support>(n);
		s3_ = new ArrayList<Support>(n);
		s4_ = new ArrayList<Support>(n);
		
		for (int i = 0; i < n; i++) {
			s1_.add(new Support(n));
			s2_.add(new Support(n));
			s3_.add(new Support(n));
			s4_.add(new Support(n));
		}
	}
	
	// ----------------------------------------------------------------------------
	
	/** Proceed on iteration of the numerical integration (step size = h). */
	public void advance(final double t, final double h, final DoubleMatrix1D dW,
			final DoubleMatrix1D dZ, final DoubleMatrix1D Xin, DoubleMatrix1D Xout) throws Exception {
		
		int n = system_.getDimension();
		double sqrt_h = Math.sqrt(h);
		double m = (double) n;
		double temp = 0;
		
	    // Create supporting vectors
	    for (int i = 0; i < n; i++)
	    {
	        for (int j = 0; j < n; j++)
	        {
	        	temp = Xin.get(j) + 1/m*F_.get(j)*h;
	        	s1_.get(i).supportX_.set(j, temp + G_.get(j, i)*sqrt_h);
	        	s2_.get(i).supportX_.set(j, temp - G_.get(j, i)*sqrt_h);
	        }

	        // "t" irrelevant (autonomous) here
	        system_.getDriftAndDiffusion(t, s1_.get(i).supportX_, s1_.get(i).supportF_, s1_.get(i).supportG_);
	        system_.getDriftAndDiffusion(t, s2_.get(i).supportX_, s2_.get(i).supportF_, s2_.get(i).supportG_);
	        
	        for (int j=0; j<n; j++) {
	        	s3_.get(i).supportX_.set(j, s1_.get(i).supportX_.get(j) + s1_.get(i).supportG_.get(j, i)*sqrt_h);
	        	s4_.get(i).supportX_.set(j, s1_.get(i).supportX_.get(j) - s1_.get(i).supportG_.get(j, i)*sqrt_h);
	        }
	        
	        // "t" irrelevant (autonomous) here
	        system_.getDriftAndDiffusion(t, s3_.get(i).supportX_, s3_.get(i).supportF_, s3_.get(i).supportG_);
	        system_.getDriftAndDiffusion(t, s4_.get(i).supportX_, s4_.get(i).supportF_, s4_.get(i).supportG_);
	    }
		
	    // For each dimension, iterate the solution by h.
	    // For clarity we'll do this line for line as per Kloeden et al. p380
	    for (int i = 0; i < n; i++) {
	    
	        // 1
	    	Xout.set(i, Xin.get(i) + G_.get(i, i)*dW.get(i));

	        // 2
	    	temp = 0;
	    	for (int j = 0; j < n; j++) {
	    		temp += (s1_.get(j).supportF_.get(i) - s2_.get(j).supportF_.get(i))*dZ.get(i);
	    	}
	    	temp *= 1.0/(2.0*sqrt_h);
	    	Xout.set(i, Xout.get(i) + temp);

	        // 3
	    	temp = 0;
	    	for (int j = 0; j < n; j++) {
	    		temp += s1_.get(j).supportF_.get(i) - 2.0*(m-2.0)/m*F_.get(i) + s2_.get(j).supportF_.get(i);
	    	}
	    	temp *= (1.0/4.0)*h;
	    	Xout.set(i, Xout.get(i) + temp);

	        // 4
	    	temp = 1.0/(4.0*sqrt_h) * (s1_.get(i).supportG_.get(i, i) - s2_.get(i).supportG_.get(i, i)) * (Math.pow(dW.get(i), 2.0) - h);
	    	Xout.set(i, Xout.get(i) + temp);

	        // 5
	    	temp = 1.0/(2.0*h) * (s1_.get(i).supportG_.get(i, i) - 2.0*G_.get(i, i) + s2_.get(i).supportG_.get(i, i)) * (dW.get(i)*h - dZ.get(i));
	    	Xout.set(i, Xout.get(i) + temp);

	        // 6
	    	temp = 1.0/(4.0*h) * (s3_.get(i).supportG_.get(i, i) - s4_.get(i).supportG_.get(i, i) - s1_.get(i).supportG_.get(i, i) + s2_.get(i).supportG_.get(i, i)) * (1.0/3.0*Math.pow(dW.get(i), 2.0) - h) * dW.get(i);
	    	Xout.set(i, Xout.get(i) + temp);
	    }
	}
	
	// ----------------------------------------------------------------------------
	
	/** Returns a description for this solver. */
	public String getDescription() {
		
		String desc = "Explicit stochastic Runge-Kutta method for Ito SDEs (Kloeden et al., 1994).\n";
		desc += "Strong order of convergence: 1.5\n";
		//desc += "Weak order of convergence ?\n";
		desc += "Only diagonal noise is handled.";
		
		return desc;
	}
	
	// ============================================================================
	// PRIVATE CLASSES
	
	/** Returns the description of this solver. */
	private class Support {
	
		/** Support vector. */
		public DoubleMatrix1D supportX_;
		/** Drift evaluated at the support vector. */
		public DoubleMatrix1D supportF_;
		/** Diffusion evaluated at the support vector. */
		public DoubleMatrix2D supportG_;
		
		public Support(int dimension) {
			supportX_ = new DenseDoubleMatrix1D(dimension);
			supportF_ = new DenseDoubleMatrix1D(dimension);
			supportG_ = new DenseDoubleMatrix2D(dimension, dimension);
		}
	}
}
