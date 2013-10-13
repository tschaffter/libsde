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

/** 
 * <p>Implements the explicit Euler-Maruyama method (strong order of convergence 0.5) to use with Ito scheme.</p>
 * 
 * <p>The iterative Euler-Maruyama method is given by the formula described in [1]:</p>
 * 
 * <p>Xn+1 = Xn + Fn*h + Gn*dWn (Euler)</p>
 * 
 * <p>where h is the integration step size. This method requires to use the Ito scheme.
 * If one wants to integrate using the Stratonovich scheme, the Euler-Heun method
 * has to be applied.</p>
 * 
 * <p>If the diffusion term is zero, then the problem is fully deterministic. In this case,
 * the Euler method achieves strong order of convergence 1. Following [1], 
 * Euler method can have strong order of convergence of 1 in case of additive noise,
 * i.e. when the diffusion term is constant and has the form:</p>
 * 
 * <p>G(t,x) := G(t)</p>
 * 
 * <p>with appropriate smoothness assumptions on both drift F and diffusion G term.</p>
 * 
 * <p>"Note however that the most widely used Euler scheme for the numeric solution of 
 * Langevin equations requires the equation to be in Ito form."<br>
 * http://en.wikipedia.org/wiki/Stratonovich_integral</p>
 * 
 * <p>If one wants to integrate using the Stratonovich scheme, the Euler-Heun method
 * can be used [1-2].</p>
 * 
 * <p>Xn+1 = Xn + Fn*h + 0.5[Gn + G(Xaux)]*dWn (Euler-Heun)<br>
 * Xaux = Xn + Gn*dWn</p>
 * 
 * <p>[1] P.E. Kloeden, E. Platen, and H. Schurz, Numerical solution of SDE through
 *     computer experiments, Springer, 1994. (pp 98, 150-153)</p>
 * <p>[2] H. Gilsing and T. Shardlow, SDELab: A package for solving stochastic differ-
 *     ential equations in MATLAB, Journal of Computational and Applied Math-
 *     ematics 205 (2007), no. 2, 1002-1018.</p>
 *     
 * @version June 28, 2009
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class EulerMaruyama extends SdeSolver {
	
	// ============================================================================
	// PUBLIC METHODS

	/** Default constructor. */
	public EulerMaruyama() {
		
		super();
	}
	
	// ----------------------------------------------------------------------------
	
	/** Proceed on iteration of the numerical integration (step size = h). */
	public void advance(final double t, final double h, final DoubleMatrix1D dW,
			final DoubleMatrix1D dZ, DoubleMatrix1D Xin, DoubleMatrix1D Xout) throws Exception {
		
		int n = system_.getDimension();
		
		for (int i = 0; i < n; i++)
			Xout.set(i, Xin.get(i) + F_.get(i)*h + G_.get(i, i)*dW.get(i));
	}

	// ----------------------------------------------------------------------------
	
	/** Returns the description of this solver. */
	public String getDescription() {
		
		String desc = "Explicit Euler-Maruyama method for Ito SDEs (Kloeden et al., 1994).\n";
		desc += "Strong order of convergence: 0.5\n";
		desc += "Weak order of convergence: 1\n";
		desc += "Only diagonal noise is handled.";
		
		return desc;
	}
}
