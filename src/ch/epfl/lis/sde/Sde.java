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

package ch.epfl.lis.sde;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * Describes a system of stochastic differential equations (SDE).
 * 
 * @version June 28, 2009
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public abstract class Sde {
	
	/** Number of SDE in the system. */
	protected int dimension_ = 1;
	
	/** Identifier for this system of SDEs. */
	protected String name_ = "";
	
	// ============================================================================
	// ABSTRACT METHODS
	
	/** Computes the drift F and diffusion G coefficients at a time t. */
	public abstract void getDriftAndDiffusion(final double t, final DoubleMatrix1D Xin, 
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor. */
	public Sde() {}
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Describes a new system of equations.
	 * @param dimension Number of differential equations in the system.
	 */
	public Sde(int dimension) {
		
		dimension_ = dimension;
	}
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public void setDimension(int dimension) { dimension_ = dimension; } 
	public int getDimension() { return dimension_; }
	
	public void setName(String name) { name_ = name; }
	public String getName() { return name_; }
}
