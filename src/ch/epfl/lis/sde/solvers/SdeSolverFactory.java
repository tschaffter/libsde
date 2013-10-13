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

import com.esotericsoftware.minlog.Log;

/**
 * Factory to instantiate easily SDE solvers.
 * 
 * @version June 28, 2009
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class SdeSolverFactory {

	public static final int EULER_MARUYAMA = 1;
	public static final int EULER_HEUN = 2;
	public static final int MILSTEIN_ITO = 3;
	public static final int MILSTEIN_STRATONOVICH = 4;
	public static final int SRK_ITO = 5;
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor. */
	public SdeSolverFactory() {}
	
	// ----------------------------------------------------------------------------
	
	/** Returns the desired solver. */
	public static SdeSolver createSolver(final int type) {
		
		SdeSolver solver = null;
		
		switch (type) {
			case EULER_MARUYAMA:
				solver = new EulerMaruyama();
				break;
			case EULER_HEUN:
				solver = new EulerHeun();
				break;
			case MILSTEIN_ITO:
				solver = new MilsteinIto();
				break;
			case MILSTEIN_STRATONOVICH:
				solver = new MilsteinStratonovich();
				break;
			case SRK_ITO:
				solver = new SRK15();
				break;
			default:
				Log.info("SdeSolverFactory", "Unrecognized solver.");
		}
		
		return solver;
	}
}
