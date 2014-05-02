/**
 * 
 */

package de.uni_jena.bio.informatik.functions;

import java.io.FileNotFoundException;
import java.util.Comparator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * SparseMatrix.
 * <h3>Usage</h3>
 * <ol>
 * <li>Generates a sparse matrix from a distance matrix {@code generateSparseMatrix} </li>
 * <li>Sorts the sparse matrix in ascending order</li>
 * </ol>
 * 
 * @author Purva Kulkarni
 *
 */

public class SparseMatrix{

	// Default Constructor
	public  SparseMatrix() {};

	/**
	 * Method to generate sparse matrix and sort it
	 * 
	 * @param float[][] distanceMatrix
	 * @param String[][] coordinateMatrix
	 * @return Set<Map.Entry<String, Float>> // returns a set coordinate and distance values
	 */	

	public Set<Map.Entry<String, Float>> generateSparseMatrix(float[][] distanceMatrix, String[][] coordinateMatrix)
			throws FileNotFoundException {

		TreeMap<String, Float> sparse = new TreeMap<String, Float>();
		ValueComparator bvc =  new ValueComparator(sparse);
		TreeMap<String,Float> sorted_map = new TreeMap<String,Float>(bvc);

		for(int row = 0; row < distanceMatrix.length;row++)
		{
			for(int column = 0; column < distanceMatrix.length; column++)
			{
				if(distanceMatrix[row][column] != 0)
				{
					String coordinateInDistanceMatrix = Integer.toString(row) +"\t"+ Integer.toString(column) + "\t" + coordinateMatrix[row][column];

					sparse.put(coordinateInDistanceMatrix, new Float(distanceMatrix[row][column]));

				}
			}

		}
		sorted_map.putAll(sparse);

		// Generate a sorted sparse set in ascending order and return it 
		Set<Map.Entry<String, Float>> set = sorted_map.entrySet();
		return set;
	}
}

/**
 * Comparator to sort sparse matrix
 * 
 */	

// Comparator to sort the sparse Map
class ValueComparator implements Comparator<String> {

	Map<String, Float> base;
	public ValueComparator(Map<String, Float> base) {
		this.base = base;
	}

	// Note: this comparator imposes orderings that are inconsistent with equals.    
	public int compare(String a, String b) {
		if (base.get(a) <= base.get(b)) {
			return -1;
		} else {
			return 1;
		} // returning 0 would merge keys
	}
}