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

import com.esotericsoftware.minlog.Log;

import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;

/**
 * Settings and functions used by libSDE (Singleton design pattern).
 * 
 * @version June 28, 2009
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class SdeSettings {

	/** Instance. */
	private static SdeSettings instance_ = null;
	
	/** Integration step size of the Wiener process (Brownian path). */
	private double dt_ = 0.01;
	/**
	 * Number of Wiener increments (multiplier_) in an integration step size (dt_).<br>
	 * <b>Important:</b> Do not touch this parameter instead you know what you do.
	 * */
	private int multiplier_ = 1;
	/** Solver ends at t = maxt_ (start at t=0). */
	private double maxt_ = 100.;
	
	/** Mersenne Twister random number generator. */
	private MersenneTwister mersenneTwister_ = null;
	/** Normal distribution N(0,1) (uses Polar Box-Muller transformation). */
	private Normal normalDistribution_ = null;
	/** Seed used to generate the Wiener process (set to -1 to initialize with java.util.Date()). */
	private int seed_ = -1;
	
	/** Get Singleton instance. */
	public static SdeSettings getInstance() {
		
		if (instance_ == null)
			instance_ = new SdeSettings();
		return instance_;
	}
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor. */
	public SdeSettings() {
		
		Log.setLogger(new SdeLogger());
	}
	
	// ----------------------------------------------------------------------------
	
	/** Initializes the random number generator with the specified seed. */
	public void initializeRNG() {
		
		if (seed_ == -1)
			mersenneTwister_ = new MersenneTwister(new java.util.Date());
		else
			mersenneTwister_ = new MersenneTwister(seed_);
		
		// // mean = 0 and std = 1
		normalDistribution_ = new Normal(0, 1, mersenneTwister_);
	}
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public void setDt(double dt) { dt_ = dt; }
	public double getDt() { return dt_; }
	
	public void setMultiplier(int multiplier) { multiplier_ = multiplier; }
	public int getMultiplier() { return multiplier_; }
	
	public void setMaxt(double T) { maxt_ = T; }
	public double getMaxt() { return maxt_; }
	
	public void setSeed(int seed) { seed_ = seed; }
	public int getSeed() { return seed_; }
	
	public Normal getNormalDistribution() { return normalDistribution_; }
}
