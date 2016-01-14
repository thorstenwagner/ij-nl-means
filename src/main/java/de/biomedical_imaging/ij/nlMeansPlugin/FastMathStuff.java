package de.biomedical_imaging.ij.nlMeansPlugin;

public class FastMathStuff {
	private static final int    BIG_ENOUGH_INT   = 16 * 1024;
    private static final double BIG_ENOUGH_ROUND = BIG_ENOUGH_INT + 0.5;
    
    
    
	public static double max(final double a, final double b) {
	        if (a > b) {
	            return a;
	        }
	        if (a < b) {
	            return b;
	        }
	        /* if either arg is NaN, return NaN */
	        if (a != b) {
	           return Double.NaN;
	        }
	        /* min(+0.0,-0.0) == -0.0 */
	        /* 0x8000000000000000L == Double.doubleToRawLongBits(-0.0d) */
	        long bits = Double.doubleToRawLongBits(a);
	        if (bits == 0x8000000000000000L) {
	            return b;
	        }
	        return a;
	    }
	
	    // http://www.java-gaming.org/index.php?topic=24194.0
	    public static int fastRound(double x) {
	       return (int) (x + BIG_ENOUGH_ROUND) - BIG_ENOUGH_INT;
	    }
	

}
