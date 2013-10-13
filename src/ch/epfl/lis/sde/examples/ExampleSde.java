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

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import ch.epfl.lis.sde.Sde;

/**
 * Simple SDE to show how to integrate a SDE using libSDE.
 * 
 * <p>dX = (-3X+1)dt + sigma*dW	with	X(0) = 1	(1)</p>
 * 
 * <p>with the expected solution:</p>
 * 
 * <p>E(X) = 2/3*exp(-3t) + 1/3						(2)</p>
 * 
 * <p>Ito and Stratonovich drift terms are equal since sigma is constant. We
 * remember the Ito-Stratonovich drift conversion:</p>
 * 
 * <p>Fu = F - 0.5*G*G'</p>
 * 
 * <p>where Fu is the drift term in Stratonovich scheme and G the diffusion term.</p>
 * 
 * @version June 28, 2009
 *
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class ExampleSde extends Sde {

	private double sigma_;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor that builds a system of 1 equations. */
	public ExampleSde() {
		
		this(1);
	}
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Constructor that builds a system of n identical equations.
	 * @param n Dimension of the system.
	 */
	public ExampleSde(int n) {
		
		super(n);
		name_ = "ExampleSde";
		sigma_ = 0.2;
	}

	// ----------------------------------------------------------------------------

	/**
	 * <p>Computes the drift coefficients F and diffusion coefficient G at a given time.
	 * Note that Ito and Stratonovich drift terms are equivalent because the diffusion
	 * term G is constant. See class description.</p>
	 */
	@Override
	public void getDriftAndDiffusion(final double t, final DoubleMatrix1D Xin,
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception {
		
		for (int i = 0; i < dimension_; i++) {
			F.set(i, -3*Xin.get(i) + 1);
			G.set(i, i, sigma_);
		}
	}
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public void setSigma(double sigma) { sigma_ = sigma; }
	public double getSigma() { return sigma_; }
}
