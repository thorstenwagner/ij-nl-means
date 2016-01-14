package de.biomedical_imaging.ij.nlMeansPlugin.test;

import static org.junit.Assert.*;
import ij.ImagePlus;

import java.net.URL;

import org.junit.Test;

import de.biomedical_imaging.ij.nlMeansPlugin.NLMeansDenoising_;

public class NLMeansDenoisingTest {

	@Test
	public void applyNonLocalMeansTest_SquareImage() {
		URL url = this.getClass().getClassLoader().getResource("square.png");
		ImagePlus input = new ImagePlus(url.getPath());
		NLMeansDenoising_ nlm = new NLMeansDenoising_();
		nlm.applyNonLocalMeans(input.getProcessor(), 20);
		byte[] nlOfInput = (byte[]) input.getProcessor().getPixels();
		url = this.getClass().getClassLoader().getResource("squareResult.png");
		ImagePlus result = new ImagePlus(url.getPath());
		byte[] expectedPixels = (byte[]) result.getProcessor().getPixels(); 
		assertArrayEquals(expectedPixels, nlOfInput);
	}
	
	@Test
	public void applyNonLocalMeansTest_HigherThanBroadImage() {
		URL url = this.getClass().getClassLoader().getResource("higher.png");
		ImagePlus input = new ImagePlus(url.getPath());
		NLMeansDenoising_ nlm = new NLMeansDenoising_();
		nlm.applyNonLocalMeans(input.getProcessor(), 20);
		byte[] nlOfInput = (byte[]) input.getProcessor().getPixels();
		url = this.getClass().getClassLoader().getResource("higherResult.png");
		ImagePlus result = new ImagePlus(url.getPath());
		byte[] expectedPixels = (byte[]) result.getProcessor().getPixels(); 
		assertArrayEquals(expectedPixels, nlOfInput);
	}
	
	@Test
	public void applyNonLocalMeansTest_BroaderThanHighImage() {
		URL url = this.getClass().getClassLoader().getResource("broader.png");
		ImagePlus input = new ImagePlus(url.getPath());
		NLMeansDenoising_ nlm = new NLMeansDenoising_();
		nlm.applyNonLocalMeans(input.getProcessor(), 15);
		byte[] nlOfInput = (byte[]) input.getProcessor().getPixels();
		url = this.getClass().getClassLoader().getResource("broaderResult.png");
		ImagePlus result = new ImagePlus(url.getPath());
		byte[] expectedPixels = (byte[]) result.getProcessor().getPixels(); 
		assertArrayEquals(expectedPixels, nlOfInput);
	}

}
