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

import java.util.Date;

import com.esotericsoftware.minlog.Log.Logger;

/**
 * Defines a custom logger for libSDE using the minlog logging framework.
 * 
 * <p>The implementation is inspired from the original Log.Logger class of
 * the library minlog.</p>
 * 
 * @version December 6, 2011
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 */
public class SdeLogger extends Logger {
	
	/** Reference time of the first log. */
	@SuppressWarnings("unused")
	private long firstLogTime_ = new Date().getTime();
	
	// ============================================================================
	// PUBLIC METHODS
	
    public void log(int level, String category, String message, Throwable ex) {
    	
    	StringBuilder builder = new StringBuilder(256);

//    	long time = new Date().getTime() - firstLogTime_;
//    	long minutes = time / (1000 * 60);
//    	long seconds = time / (1000) % 60;
//    	if (minutes <= 9) builder.append('0');
//		builder.append(minutes);
//		builder.append(':');
//		if (seconds <= 9) builder.append('0');
//		builder.append(seconds);

//		switch (level) {
//			case Log.LEVEL_ERROR:
//				builder.append(" ERROR: ");
//				break;
//			case Log.LEVEL_WARN:
//				builder.append("  WARN: ");
//				break;
//			case Log.LEVEL_INFO:
//				builder.append("  INFO: ");
//				break;
//			case Log.LEVEL_DEBUG:
//				builder.append(" DEBUG: ");
//				break;
//			case Log.LEVEL_TRACE:
//				builder.append(" TRACE: ");
//				break;
//		}

		if (category != null) {
			builder.append('[');
			builder.append(category);
			builder.append("] ");
		}
		
		builder.append(message);
        
        // print to the standard output
        System.out.println(builder);
        
        if (ex != null)
        	ex.printStackTrace();
    }
}