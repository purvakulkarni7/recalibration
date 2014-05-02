/**
 * 
 */

package de.uni_jena.bio.informatik.functions;

/**
 * MatrixFunctions.
 * <h3>Usage</h3>
 * <ol>
 * <li>Transposes a two dimensional float matrix {@code transposeFloat} </li>
 * <li>Transposes a two dimensional float matrix {@code transposeString}</li>
 * <li>Find the minimum value in a two dimensional float matrix and 
 * and returns the coordinate position of this value {@code minValue} </li>
 * </ol>
 * 
 * @author Purva Kulkarni
 *
 */

public class MatrixFunctions {

	/**
	 * Method to transpose double matrix
	 * 
	 * @param float[][] matrix
	 * @param int size
	 * @return float[][] transposeDouble
	 */	
	
	public float[][] transposeFloat(float[][] matrix, int size)
	{
		int c, d;
		float transposeFloat[][] = new float[size][size];

		for ( c = 0 ; c < size ; c++ )
		{
			for ( d = 0 ; d < size ; d++ )               
				transposeFloat[d][c] = matrix[c][d];
		}


		return transposeFloat;
	}

	/**
	 * Method to transpose string matrix
	 * 
	 * @param String[][] matrix
	 * @param int size
	 * @return String[][] transposeString 
	 */	

	public String[][] transposeString(String[][] matrix, int size)
	{
		int c, d;
		String transposeString[][] = new String[size][size];

		for ( c = 0 ; c < size ; c++ )
		{
			for ( d = 0 ; d < size ; d++ )               
				transposeString[d][c] = matrix[c][d];
		}


		return transposeString;
	}

	/**
	 * Method to calculate minimum value and its coordinate position in the matrix
	 * 
	 * @param float[][] matrix
	 * @return int coordinateI 
	 */	
	
	float minValue = 0;
	public int minValue(float[][] matrix) 
	{

		for (int i = 0; i < matrix.length; i++) 
		{
			for (int j = 0; j < matrix[i].length; j++) 
			{
				if(matrix[i][j] != 0 )
				{
					minValue = matrix[i][j]; // Assign the first non-zero value as minValue
					break;
				}
				else
					continue;
			}
		}
		
		int i,j;
		int coordinateI = 0, coordinateJ = 0 ;
		for (i = 0; i < matrix.length; i++) 
		{
			for (j = 0; j < matrix[i].length; j++) 
			{
				if (matrix[i][j] < minValue && matrix[i][j] !=0) 
				{
					minValue = matrix[i][j];
					coordinateI = i;
					coordinateJ = j;

				}
			}
		}
		System.out.println("Min Value is: " + minValue + "Coordinates: " + coordinateI + " " + coordinateJ );

		return coordinateI;
	}
}
