/**
 * 
 */

package de.uni_jena.bio.informatik.main;

import java.io.IOException;

/**
 * SpectralRecalibration.
 * <h3>Usage</h3>
 * <ol>
 * <li>Class containing main method.</li>
 * <li>Creates an InputData object {@code InputData} and calls the userInput {@code userInput} 
 * to ask for user input</li>
 * </ol>
 * 
 * @author Purva Kulkarni
 *
 */

import de.uni_jena.bio.informatik.input.InputData;

public class SpectralRecalibration 
{
	// main method begins execution of java application
	public static void main(String[] args) throws IOException {
		//InputData f1 = new InputData();
		InputData f1 = new InputData();

		f1.userInput();

	} //end main method

} // end class
