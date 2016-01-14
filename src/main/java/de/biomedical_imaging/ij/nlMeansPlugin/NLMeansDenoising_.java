/**
 * Non local means filter for ImageJ
 * Copyright (C) 2013  Pascal Behnel & Thorsten Wagner
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
package de.biomedical_imaging.ij.nlMeansPlugin;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Pascal
 */
public class NLMeansDenoising_ implements PlugInFilter {
    
    private final int weightPrecision = 1000; // Precision of the used Weight. Values > 1000 can cause Overflow

    private ImagePlus imp; // ImagePlus representation of the Image
    private int width; // Width of the Image
    private int height; // Height of the Image
    
    private int[][] pixels; // Pixels of the Image. First Dimension is for Colour-Channel, second for Pixels. Pixels are arranged in width * height
//    private int[][] pixelsExpand; // Pixels of the expanded Version of the Image. This Version is needed to prevent Walking out of Bounds.
    private int widthE; // Width of the expanded Image
    private int heightE; // Height of the expanded Image
    
    private int w; // Big Search-Window
    private int n; // Small Search-Window used for Patches
    
    private double sigma2; // Variance of the Image-Noise
    private double h2; // Smoothing-Parameter
    private int distConst; // Constant Value that the Distance gets Multiplied with. Currently unused.
    private int dim; // Dimension of the Image. (1 = Grayscale, 3 = RGB)
    
    private int nextdx; // Variable used to store the next Value for dx. 
    private int nextdy; // Variable used to store the next Value for dy. 
    
    private long[][] uL; // Long Representation of the Denoising Image
    private long[] wMaxArrL; // Max-Weight that was used for each Pixel. (Long Representation)
    private long[] wSumArrL;// Sum of all Weights per Pixel (Long Representation)
    
    private boolean autoEstimate = false; // Use Auto-Estimation to determine Image-Noise
    private int constantSigma = 15; // Standard-Value for Sigma
    private int smoothingFactor = 1;
    private int[] usedSigmas; //Saves the used sigmas when processing stacks..
    
    @Override
    public int setup(String arg, ImagePlus imp) {
    	
    	if (arg.equals("final")) {
    		//Save used sigmas
    		String sigmas = "" + usedSigmas[0];
    		for(int i = 1; i < usedSigmas.length; i++){
    			sigmas += "," + usedSigmas[i];
    		}
    		Prefs.set("nlmeans.sigma", sigmas);
    		Prefs.set("nlmeans.smoothingfactor", smoothingFactor);
			return DONE;
		}
    	
        this.imp = imp;
//        this.ip = imp.getProcessor();
        GenericDialog gd = new GenericDialog("Non-Local Means");
        gd.addNumericField("Sigma", 15, 0);
//        gd.addNumericField("Window Width", 512, 0);
//        gd.addNumericField("Window Height", 512, 0);
        gd.addNumericField("Smoothing_Factor", 1, 0);
        gd.addCheckbox("Auto estimate sigma", false);
        gd.addHelp("http://fiji.sc/Non_Local_Means_Denoise");
        gd.showDialog();
        if (gd.wasCanceled()) {
            return -1;
        }
        constantSigma = (int)gd.getNextNumber();
        smoothingFactor = (int)gd.getNextNumber();
        autoEstimate = gd.getNextBoolean();
        usedSigmas = new int[imp.getNSlices()*imp.getNFrames()];
        
        return IJ.setupDialog(imp, DOES_ALL+FINAL_PROCESSING);
    }
    
    @Override
    public void run(ImageProcessor ip) {
        int sigma = constantSigma;
        if (autoEstimate) {
            sigma = (int)getGlobalNoiseLevel(ip);
        } 
        usedSigmas[ip.getSliceNumber()-1] =sigma; 
        sigma = smoothingFactor*sigma;
        
        applyNonLocalMeans(ip, sigma);
   
    }
    
    /**
     * 
     * @param ip Image which should be denoised
     * @param sigma Estimated standard deviation of noise
     * @param imageType Type of image (e.g. ImagePlus.Gray8 etc.)
     */
    public void applyNonLocalMeans(ImageProcessor ip, int sigma){
        initSettings(sigma, ip);
        
        try {
          int width = 512;
          int height = 512;
          double[][] result = NLMeansDenoising(ip, width, height);
          createPicture(result, ip);
        } catch (InterruptedException e) {
    	  e.printStackTrace();
      //    IJ.showMessage("Error while computing Denoised Image.");
      }
    }

    
    private double[][] NLMeansDenoising(ImageProcessor ip, int windowWidth, 
            int windowHeight) throws InterruptedException {
        
        double[][] result = new double[dim][ip.getWidth()*ip.getHeight()];
        
        for (int ys = 0; ys < ip.getHeight(); ys+=windowHeight) {
            for (int xs = 0; xs < ip.getWidth(); xs+=windowWidth) {
                int imagePartWidth = (windowWidth + xs > ip.getWidth()) 
                        ? windowWidth - ((windowWidth + xs) - ip.getWidth()) : windowWidth;
                int imagePartHeight = (windowHeight + ys > ip.getHeight()) 
                        ? windowHeight - ((windowHeight + ys) - ip.getHeight()) : windowHeight;
                int[][] imagePartE = expandImage(pixels, xs, ys, 
                        imagePartWidth, 
                        imagePartHeight, 
                        ip.getWidth(), ip.getHeight(), false);
                
//                double[][] partResult = NLMeansMultithreadInstance(imagePartE, 
//                        Runtime.getRuntime().availableProcessors(), imagePartWidth, 
//                        imagePartHeight);
                double[][] partResult = NLMeansMultithreadInstance(imagePartE, 
                        1, imagePartWidth, 
                        imagePartHeight);
                
                // save Partial Result in Image
                nextdx = -w; 
                nextdy = -w;
                int ystart = ys;
                int xstart = xs;
                for (int y = ystart; y < ystart + imagePartHeight; y++) {
//                    if (y >= ip.getHeight()) continue;
                    for (int x = xstart; x < xstart + imagePartWidth; x++) {
//                        if (x >= ip.getWidth()) continue;
                        for (int d = 0; d < dim; d++) {
                            result[d][y*ip.getWidth() + x] = partResult[d][(y - ystart)*imagePartWidth + x - xstart];
                        }
                    }
                }
            }
        }
        
        
        
        
        
        return result;
    }    

    /**
     * Multi Threaded Implementation of the Non-local Means Algorithm.
     * This accelerated Version is based of: Darbon, Jérôme, et al. 
     * "Fast nonlocal filtering applied to electron cryomicroscopy." 
     * Biomedical Imaging: From Nano to Macro, 2008. ISBI 2008. 
     * 5th IEEE International Symposium on. IEEE, 2008.
     * @param image The image as Integer Array. Colors are stored within first 
     * dimension of Array. Gets computed via convertImage()
     * @param threadcount Number of Threads used for Denoising
     * @param ip ImageProcessor for the original Image
     * @throws InterruptedException 
     */
    private double[][] NLMeansMultithreadInstance(int[][] image, int threadcount,
            int width, int height) throws InterruptedException {
        int widthE = width + 2*w + 2*n;
        int heightE = height + 2*w + 2*n;
        long[][] u = new long[dim][widthE * heightE];
        long[] wMaxArr = new long[widthE * heightE];
        long[] wSumArr = new long[widthE * heightE];
        
        List<Worker> workerList = new ArrayList<Worker>(threadcount);
        for (int i = 0; i < threadcount; i++) {
            Worker worker = new Worker(width, height, image, u, wMaxArr, wSumArr);
            worker.start();
            workerList.add(worker);
        }
        for (Worker worker : workerList) {
            worker.join();
        }
        
        return finishPicture(u, image, wMaxArr, wSumArr, width, height);
    }
    
    private synchronized void deliverImagePart(long[][] imagePart, long[][] u, int widthE, int heightE) {
        for (int y = 0; y < heightE; y++) {
            int offset = y*widthE;
            for (int x = 0; x < widthE; x++) {
                for (int d = 0; d < dim; d++) {
                    u[d][offset + x] += imagePart[d][offset + x];
                }
            }  
        }
    }
    
    /**
     * This Method is used to deliver a partial result of the Weight Sum Array.
     * The Weight Sum Array stores the sum of all Weights that are used
     * for each pixel. It is used within finishPicture(...) to properly Scale
     * each Pixel.
     * @param arr Weight Sum Array
     */
    private synchronized void deliverWSumArr(long[] arr, long[] wSumArr, int widthE, int heightE) {
        for (int y = 0; y < heightE; y++) {
            int offset = y*widthE;
            for (int x = 0; x < widthE; x++) {
                wSumArr[offset + x] += arr[offset + x];
            }  
        }
    } 
    
    /**
     * This Method is used to deliver a partial result of the Weight Max Array.
     * The Weight Max Array stores the maximum Weight that is used per Pixel.
     * This Weight is used as Weight between the Pixel and itself.
     * @param arr Maximum Weight Array
     */
    private synchronized void deliverWMaxArr(long[] arr, long[] wMaxArr, int widthE, int heightE) {
        for (int y = 0; y < heightE; y++) {
            int offset = y*widthE;
            for (int x = 0; x < widthE; x++) {
                if (wMaxArr[offset + x] < arr[offset + x]) {
                    wMaxArr[offset + x] = arr[offset + x];
                }
            }  
        }
    }
    
    /**
     * Finishes the Picture by dividing every Pixel with the Sum of all Weights
     * for the respective Pixel, and by performing the last denoising step.
     * As last Step, the Pixels get weighted with the maximum Weight for each Pixel.
     * @param picture The Denoised Picture
     * @param wMaxArr Array with highest used Weight for each Pixel
     * @param wSumArr Array with Sum of Weights for each Pixel
     * @return 
     */
    private double[][] finishPicture(long[][] picture, int[][] pixelsExpand, 
            long[] wMaxArr, long[] wSumArr, int width, int height) {
        double[][] result = new double[dim][width*height];
        int wn = w + n;
        int widthE = width + 2*wn;

        // x and y coordinates are based off the original Image (NOT the expanded Image)
        for (int y = 0; y < height; y++) {
            int offset = y*width; // y offset for original Image coordinates
            int offset2 = (y + wn) * widthE; // y offset for expanded Image coordinates
            for (int x = 0; x < width; x++) {
                int k = offset + x; // array Position for Pixel x, y
                int kwn = offset2 + x + wn; // same as k, but in expanded Image coordinates
                for (int d = 0; d < result.length; d++) {
                    result[d][k] = picture[d][kwn];
                    
                    if (wMaxArr[kwn] == 0) {
                        // If Sum of all Weights is 0, just copy the original Pixel
                        result[d][k] += pixelsExpand[d][kwn];
                    } else {
                        // Weight the original Pixel with the maximum Weight
                        result[d][k] += pixelsExpand[d][kwn] * wMaxArr[kwn];
                        wSumArr[kwn] += wMaxArr[kwn];
                        
                        // Divide Pixel by sum of all Weights
                        result[d][k] /= wSumArr[kwn];
                    }
                }
            }
        }
        
        return result;
    }
    
    private void denoise(long[][] targetArr, int[][] pixelsExpand, long[][] S, 
            long[] wMaxArr, long[] wSumArr, int widthE, int heightE, int dx, int dy) {
        int wn = w + n;
        for (int y = wn; y < heightE - wn; y++) {
            int offset = y*widthE;
            int offsetn = (y + dy)*widthE;
            for (int x = wn; x < widthE - wn; x++) {
                int k = offset + x;
                int kn = offsetn + x + dx;
                int weight = computeWeight(S, widthE, x, y, weightPrecision);
                wMaxArr[k] = Math.max(weight, wMaxArr[k]);
                wSumArr[k] += weight;
                wMaxArr[kn] = Math.max(weight, wMaxArr[kn]);
                wSumArr[kn] += weight;
                
                for (int d = 0; d < dim; d++) {
                    int wk = weight * pixelsExpand[d][k];
                    int wkn = weight * pixelsExpand[d][kn];

                    targetArr[d][k] += wkn;
                    targetArr[d][kn] += wk;
                }
            }
        }
    }
    
    /**
     * Computes the Weight between the Pixel x,y and the Pixel that lies
     * at x + dx, y + dy. dx and dy are implicitly given because the
     * Difference Image is based on them.
     * @param S Difference Image for a dx / dy pair
     * @param x x-Coordinate of the current Pixel
     * @param y y-Coordinate of the current Pixel
     * @param precision Precision of the Weight. Should be multiple of 10
     * @return 
     */
    private int computeWeight(long[][] S, int widthE, int x, int y, int precision) {
        double distance = computeDistance(S, widthE, x, y);
       
        double exp = Math.max(distance-sigma2, 0.0);
        
//        exp /= h2;
//        double weight = Math.exp(-exp);
        double weight = h2 / (h2 + exp);
        
//        int iWeight = FastMathStuff.fastRound(weight * precision) + 1;
//        if (iWeight == 0) iWeight = 1;

        return FastMathStuff.fastRound(weight * precision);
    }
    
    /**
     * Computes the Difference between the Surroundings of the Pixel x,y and the 
     * Pixel that lies at x + dx, y + dy. dx and dy are implicitly given 
     * because the Difference Image is based on them.
     * Is used to compute the Weights. 
     * @param S Difference Image for a dx / dy pair
     * @param x x-Coordinate of the current Pixel
     * @param y y-Coordinate of the current Pixel
     * @return 
     */
    private double computeDistance(long[][] S, int widthE, int x, int y) {
        double distance = 0;
        for (int d = 0; d < dim; d++) {
            distance += S[d][(y + n) * widthE + (x + n)] 
                      + S[d][(y - n) * widthE + (x - n)]
                      - S[d][(y - n) * widthE + (x + n)]
                      - S[d][(y + n) * widthE + (x - n)];
        }
        
        return distance;
    }
    
    /**
     * Computes the Difference Image for a given dx / dy Pair. As dx and dy can
     * be negative, image needs to be expanded to prevent out of bounds errors.
     * @param image Expanded Version of Original Image
     * @param targetArr Target Array in which the Difference Image gets stored into
     * @param dx
     * @param dy 
     */
    private void computeDifferenceImage(int[][] image, long[][] targetArr, 
            int dx, int dy, int widthE, int heightE) {
        int wn = w + n;
        long temp;
        
        // Compute very first Pixel of Image (x = 0; y = 0)
        for (int d = 0; d < dim; d++) {
            temp = image[d][wn * widthE + wn]
                        - image[d][(wn + dy)*widthE + dx + wn];
            targetArr[d][wn * widthE + wn] = temp * temp;
        }
        
        // Compute first Row of Image (y = 0)
        int offset = wn * widthE;
        int offsetdy = (wn + dy) * widthE;
        for (int x = wn + 1; x < widthE; x++) {
            for (int d = 0; d < dim; d++) {
                temp = image[d][offset + x] - image[d][offsetdy + x + dx];
                targetArr[d][offset + x] = targetArr[d][offset + x - 1] 
                        + temp * temp;
            }
        }
        
        // Compute first Column of Image (x = 0)
        for (int y = wn + 1; y < heightE; y++) {
            int offsety = y*widthE;
            offsetdy = (y + dy)*widthE;
            for (int d = 0; d < dim; d++) {
                temp = image[d][offsety + wn] - image[d][offsetdy + wn + dx];
                targetArr[d][offsety + wn] = targetArr[d][offsety - widthE + wn] + temp * temp;
            }
        }
        
        // Compute rest of the Image
        for (int y = wn + 1; y < heightE; y++) {
            offset = y*widthE;
            int offset2 = (y + dy)*widthE;
            for (int x = wn + 1; x < widthE; x++) {
                for (int d = 0; d < dim; d++) {
                    targetArr[d][offset + x]  = targetArr[d][offset + x - 1];
                    targetArr[d][offset + x] += targetArr[d][offset + x - widthE];
                    targetArr[d][offset + x] -= targetArr[d][offset + x - 1 - widthE];
                    
                    temp = image[d][offset + x] - image[d][offset2 + x + dx];
                    double temp2 = temp * temp;
                    targetArr[d][offset + x] += temp2;
                }
            }
        }
    }

    
    /**
     * Expands the boundaries of an image in all four directions. The new content
     * of the Image gets filled with the adjacent parts of the Image.
     * To view a Preview of this Image, use display = true
     * @param image Original Image
     * @param display Display Preview of generated Image
     * @return 
     */
    private int[][] expandImage(int[][] image, int xstart, int ystart, int width, int height, 
            int orgWidth, int orgHeight, boolean display) {
        int heightE = height + 2*w + 2*n;
        int widthE = width + 2*w + 2*n;
        int[][] result = new int[dim][widthE * heightE];

        for (int y = 0; y < heightE; y++) {
            int yr = y - w - n + ystart;
            
//            if (yr >= orgHeight) yr = (ystart - w - n) + yr - orgHeight;
            if (yr >= orgHeight) yr = yr - orgHeight;
            if (yr < 0) yr = height + yr;
            
            
            int offset = y * widthE;
            int offsetr = yr * orgWidth;
            for (int x = 0; x < widthE; x++) {
                int xr = x + (xstart - w - n);
//                if (xr >= orgWidth) xr = xstart + xr - orgWidth;
                if (xr >= orgWidth) xr = xr - orgWidth;
                if (xr < 0) xr = width + xr;
                for (int d = 0; d < dim; d++) {
                    result[d][offset + x] = image[d][offsetr + xr];
                }
            }
        }
        
        if (display) {
            int[] pixelsPicture = new int[result[0].length];
        
            for (int y = 0; y < heightE; y++) {
                int offset = y*widthE;
                for (int x = 0; x < widthE; x++) {
                    int p = offset + x;
                    int red = (int)result[0][p];
                    int green = (int)result[1][p];
                    int blue = (int)result[2][p];
                    int pixel = ((red & 0xff)<<16) 
                              + ((green & 0xff)<<8) 
                              + (blue & 0xff);
                    pixelsPicture[p] = pixel;
                }
            }
            
            BufferedImage bimg = convertToImage(widthE, heightE, pixelsPicture);
            ImagePlus imp2 = new ImagePlus("Expanded Image", bimg);
            imp2.show();
        }
        
        return result;
    }
    
    /**
      * Implements the gaussian noise level estimation algorithm of 
      * Immerkaer, J., 1996. Fast noise variance estimation. 
      * Computer Vision and Image Understanding, 
      * 64(2), pp.300–302. 
      * @param imp
      * @return noise level
      */
    public static double getGlobalNoiseLevel(ImageProcessor ip){
        FloatProcessor fp = null;
        
        switch(ip.getBitDepth()) {
            
            case 8:
                ByteProcessor bp = (ByteProcessor)ip;
                fp = bp.duplicate().convertToFloatProcessor();
                break;
            case 24:
                ColorProcessor cp = (ColorProcessor)ip;
                fp = cp.duplicate().convertToFloatProcessor();
                break;
            case 16:
                ShortProcessor sp = (ShortProcessor)ip;
                fp = sp.duplicate().convertToFloatProcessor();
                break;
            case 32:
                fp = (FloatProcessor) ip.duplicate();
                break;
            default:
                break; 
        }

        Convolver convolver = new Convolver();
        float[] kernel = {1,-2,1,-2,4,-2,1,-2,1};
        convolver.convolve(fp, kernel, 3, 3);

        int w = fp.getWidth();
        int h = fp.getHeight();
        double sum = 0;

        for(int x = 1; x < (w-1); x++){
            for(int y = 1; y < (h-1); y++){
                sum += Math.abs(fp.getPixelValue(x, y));
            }
        }
        double sigma = Math.sqrt(Math.PI/2)*1.0/(6.0*(w-2)*(h-2))*sum;

        return sigma;
    }
    
    
    /**
     * Initialize needed Settings
     * @param sigma An estimate of the standard deviation of the Noise-Level 
     * within the Image
     * @param ip The Image-Processor of the original Image
     */
    private void initSettings(int sigma, ImageProcessor ip) {
        int type = new ImagePlus(null, ip).getType();
        
        // Init recommended Algorithm Settings
        double hfactor;
        if (type == ImagePlus.COLOR_256 || type == ImagePlus.COLOR_RGB) {
            
            // Color Image
            
            if (sigma > 0 && sigma <= 25) {
                n = 1;
                w = 10;
//                n = 3;
//                w = 17;
                hfactor = 0.55;
            } else if (sigma > 25 && sigma <= 55) {
                n = 2;
                w = 17;
                hfactor = 0.4;
            } else {
                n = 3;
                w = 17;
                hfactor = 0.35;
            }
        } else {
            
            // Gray Image
            
            if (sigma > 0 && sigma <= 15) {
                n = 1;
                w = 10;
                hfactor = 0.4;
            } else if (sigma > 15 && sigma <= 30) {
                n = 2;
                w = 10;
                hfactor = 0.4;
            } else if (sigma > 30 && sigma <= 45) {
                n = 3;
                w = 17;
                hfactor = 0.35;
            } else if (sigma > 45 && sigma <= 75) {
                n = 4;
                w = 17;
                hfactor = 0.35;
            } else {
                n = 5;
                w = 17;
                hfactor = 0.3;
            }
        }
        
        width = ip.getWidth();
        height = ip.getHeight();
        widthE = width + 2*w + 2*n;
        heightE = height + 2*w + 2*n;        

        // Read Pixels from ImageProcessor and store them in pixels Array
        convertPixels(ip, type);

        double h = hfactor * sigma;
        sigma2 = sigma * sigma * 2 * (dim * (2 * n + 1) * (2 * n + 1));
//        sigma2 = 2 * sigma * sigma;
        distConst = (dim * (2 * n + 1) * (2 * n + 1));
//        h2 = (h * h) / (dim * (2 * n + 1) * (2 * n + 1));
        h2 = (h * h);
        
        // Multithreadding related Initializations
        nextdx = -w; 
        nextdy = -w;
    }
    
    /**
     * Returns next dx / dy Pair
     * dx and dy are needed to compute a specific iteration of the Algorithm.
     * This method provides the next unused dx / dy Pair to be used in a 
     * denoising Thread.
     * @return dx and dy as int array, in this respective order
     */
    private synchronized int[] getNextDV() {
        if (nextdy > 0) return null;
        
        int[] result = new int[] { nextdx, nextdy };

        if (nextdx == w) {
            nextdy++;
            nextdx = -w;
        } else {
            nextdx++;
        }
        
        return result;
    }
    
    /**
     * Converts the Image into its proper form sothat it can be used by the 
     * Algorithm
     * @param ip
     * @param type Type of the Image based on ImageJ ImageTypes
     */
    private void convertPixels(ImageProcessor ip, int type) {
        switch(type) {
            
            case ImagePlus.COLOR_256:
                convertColor256(ip);
                break;
            case ImagePlus.COLOR_RGB:
                convertRGB(ip);
                break;
            case ImagePlus.GRAY16:
                convertGray16(ip);
                break;
            case ImagePlus.GRAY32:
                convertGray32(ip);
                break;
            case ImagePlus.GRAY8:
                convertGray8(ip);
                break;
            default:
                break; 
        }
    }
    
    private void convertColor256(ImageProcessor ip) {
        dim = 1;
        
        byte[] pixelArray = (byte[]) ip.getPixels();
        pixels = new int[dim][width*height];
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                pixels[0][pos] = pixelArray[pos] & (0xff);
            }
        }
    }
    
    private void convertRGB(ImageProcessor ip) {
        dim = 3;
        
        int[] pixelArray = (int[])ip.getPixels();
        pixels = new int[dim][width*height];
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int qtemp = pixelArray[offset + x];
                pixels[0][offset + x] = ((qtemp & 0xff0000)>>16);
                pixels[1][offset + x] = ((qtemp & 0x00ff00)>>8);
                pixels[2][offset + x] = ((qtemp & 0x0000ff));
            }
        }
    }
    
    private void convertGray32(ImageProcessor ip) {
        dim = 1;
        
        float[] pixelArray = (float[]) ip.getPixels();
        pixels = new int[dim][width*height];
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                pixels[0][pos] = (int)pixelArray[pos];
            }
        }
    }
    
    private void convertGray16(ImageProcessor ip) {
        dim = 1;
        
        short[] pixelArray = (short[]) ip.getPixels();
        pixels = new int[dim][width*height];
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                pixels[0][pos] = (int)(pixelArray[pos] & (0xffff));
            }
        }
    }
    
    private void convertGray8(ImageProcessor ip) {
        dim = 1;
        
        byte[] pixelArray = (byte[]) ip.getPixels();
        pixels = new int[dim][width*height];
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                pixels[0][pos] = (int)(pixelArray[pos] & (0xff));
            }
        }
    }
    
    /**
     * Converts a denoised Picture back to its original Format and saves it
     * in the ImageProcessor
     * @param image
     * @param ip 
     */
    private void createPicture(double[][] image, ImageProcessor ip) {
    	
        switch(ip.getBitDepth()) {
            
            case 8:
            	createPicture8Bit(image, ip);
                break;
            case 24:
                createPicture24Bit(image, ip);
                break;
            case 16:
                createPicture16Bit(image, ip);
                break;
            case 32:
                createPicture32Bit(image, ip);
                break;
            default:
                break; 
        }
        
//        imp.repaintWindow();

//        impNew.setTitle(newImageTitle);
//        impNew.show();
        
    }
    
  
    
    private void createPicture24Bit(double[][] image, ImageProcessor ip) {
//        impNew = imp.createImagePlus();
//        impNew.setProcessor(imp.getProcessor().duplicate());
        int[] pixelsPicture = (int[])ip.getPixels();
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int p = offset + x;
                int red = (int)image[0][p];
                int green = (int)image[1][p];
                int blue = (int)image[2][p];
                int pixel = ((red & 0xff)<<16) 
                          + ((green & 0xff)<<8) 
                          + (blue & 0xff);
                pixelsPicture[p] = pixel;
            }
        }
        ip.setPixels(pixelsPicture);
    }
    
    private void createPicture32Bit(double[][] image, ImageProcessor ip) {
//        impNew = imp.createImagePlus();
//        impNew.setProcessor(imp.getProcessor().duplicate());
        float[] pixelsPicture = (float[])ip.getPixels();
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                float pixel = (float)(image[0][pos]);
                pixelsPicture[pos] = pixel;
            }
        }
        
        ip.setPixels(pixelsPicture);
    }
    
    private void createPicture16Bit(double[][] image, ImageProcessor ip) {
//        impNew = imp.createImagePlus();
//        impNew.setProcessor(imp.getProcessor().duplicate());
        short[] pixelsPicture = (short[])ip.getPixels();
        
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                short pixel = (short)(image[0][pos]);
                pixelsPicture[pos] = pixel;
            }
        }
        
        ip.setPixels(pixelsPicture);
    }
    
    private void createPicture8Bit(double[][] image, ImageProcessor ip) {
//        ImagePlus impNew = imp.createImagePlus();
//        impNew.setProcessor(imp.getProcessor().duplicate());
        byte[] pixelsPicture = (byte[])ip.getPixels();
        for (int y = 0; y < height; y++) {
            int offset = y*width;
            for (int x = 0; x < width; x++) {
                int pos = offset + x;
                byte pixel = (byte)(image[0][pos]);
                pixelsPicture[pos] = pixel;
            }
        }
        
        ip.setPixels(pixelsPicture);
    }
    
    public static BufferedImage convertToImage(int width, int height, int[] pixels) {
        int wh = width * height;
        int[] newPixels = new int[wh*3];
        for (int i = 0; i < wh; i++) {
            int rgb = pixels[i];
            int red = (rgb >> 16) & 0xFF;
            int green = (rgb >> 8) & 0xFF;
            int blue = rgb & 0xFF;
            newPixels[i*3] = red;
            newPixels[i*3 + 1] = green;
            newPixels[i*3 + 2] = blue;
        }
        
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        WritableRaster raster = (WritableRaster) image.getData();
        raster.setPixels(0,0,width,height,newPixels);
        image.setData(raster);
        
        return image;
    }
    
    class Worker extends Thread {
        private int[][] image;
        
        private long[][] u;
        private long[] wMaxArr;
        private long[] wSumArr;
        
        int width;
        int height;
        
        public Worker() {
        }
        
        public Worker(int width, int height, int[][] image, long[][] u, long[] wMaxArr, long[] wSumArr) {
            this.width = width;
            this.height = height;
            this.image = image;
            this.u = u;
            this.wMaxArr = wMaxArr;
            this.wSumArr = wSumArr;
        }

        @Override
        public void run() {
            int[] vec;
            int dx, dy;
            int heightE = height + 2*w + 2*n;
            int widthE = width + 2*w + 2*n;
            long[] TwMaxArr = new long[widthE*heightE];
            long[] TwSumArr = new long[widthE*heightE];
            long[][] TimagePart = new long[dim][widthE*heightE];
            long[][] TS = new long[dim][widthE*heightE];

            vec = getNextDV();
            while (vec != null) {
                dx = vec[0];
                dy = vec[1];
                if ((2*w + 1) * dy + dx >= 0) {
                    vec = getNextDV();
                    continue;
                }
                
                // compute Sdx
                computeDifferenceImage(image, TS, dx, dy, widthE, heightE);

                // denoise with Sdx
                denoise(TimagePart, image, TS, TwMaxArr, TwSumArr, widthE, heightE, dx, dy);
                
                // get next Vector
                vec = getNextDV();
            }
            
            // save to global variables
            deliverImagePart(TimagePart, u, widthE, heightE);
            deliverWMaxArr(TwMaxArr, wMaxArr, widthE, heightE);
            deliverWSumArr(TwSumArr, wSumArr, widthE, heightE);
        }
    }
}
