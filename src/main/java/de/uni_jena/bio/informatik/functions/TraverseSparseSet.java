/**
 * 
 */

package de.uni_jena.bio.informatik.functions;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.analysis.UnivariateFunction;
import de.uni_jena.bio.informatik.input.InputData;
import de.uni_jena.bio.informatik.input.RawDataStorage;
import de.unijena.bioinf.ChemistryBase.ms.MutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import de.unijena.bioinf.recal.MzRecalibration;

/**
 * TraverseSparseSet.
 * <h3>Usage</h3>
 * <ol>
 * <li>Traverse the sparse set containing coordinate pairs and distance values</li>
 * <li>Takes the two mass lists with minimum distance and performs recalibration 
 * using functions available in the recalibration module {@code maxLinePairStabbing}</li>
 * </ol>
 * 
 * @author Purva Kulkarni
 *
 */

public class TraverseSparseSet{
	SimpleSpectrum newConsensus;
	int counter;
	// Default Constructor
	public  TraverseSparseSet() {};
	public InputData id = new InputData();

	/**
	 * Method to traverse the sparse set and perform recalibration using the recalibration module 
	 * 
	 * @param Set<Map.Entry<String, Float>> sortedSet
	 */	

	public SimpleSpectrum traverseSortedSet(Set<Map.Entry<String, Float>> sortedSet)
	{	
		TreeMap<Integer, String> rankList = generateRankList(sortedSet);

		for(int i = 0; i < rankList.size(); i++)
		{
			if(i == 0)
			{
				newConsensus = generateConsensus(rankList.get(i+2), rankList.get(i+1));
				i++;
			}
			else
			{
				System.out.println("value of i: " + i + " and taking rank :" + (i+1));
				newConsensus = generateConsensus(newConsensus, rankList.get(i+1));
			}

		}
		return newConsensus;
	}



	/** Generates a RankList when a sorted Set is provided 
	 * @param sortedSet
	 * @return TreeMap with the Rank and the file name(Rank= *.txt)
	 */
	public TreeMap<Integer, String> generateRankList(
			Set<Map.Entry<String, Float>> sortedSet) {
		String recalibrateJFileName;
		String consensusIFileName;
		int rank =0;
		TreeMap<Integer, String> rankList = new TreeMap<Integer, String>();
		for(Map.Entry<String, Float> element : sortedSet) 
		{
			String Coordinate = element.getKey();
			String delimiter = "\\t";
			String var[] = Coordinate.split(delimiter);

			//	String consensusIFile = var[0];		
			//	String recalibrateJFile = var[1];			
			//	String recalibrateJCoordinateX = var[2];				
			//	String recalibrateJCoordinateY = var[4];		
			recalibrateJFileName = var[5];		
			//	String consensusICoordinateX = var[6];		
			//	String consensusICoordinateY = var[8];
			consensusIFileName = var[9];
			//if(spectraTraversed.isEmpty())
			if(rankList.isEmpty())
			{
				rank++;
				rankList.put(rank, consensusIFileName);
				rank++;
				rankList.put(rank,recalibrateJFileName);
			}
			else
			{

				if((rankList.containsValue(consensusIFileName)) && (!rankList.containsValue(recalibrateJFileName)))
				{	
					rank++;
					rankList.put(rank,recalibrateJFileName);
				}
				else
					if((!rankList.containsValue(consensusIFileName)) && (rankList.containsValue(recalibrateJFileName)))
					{
						rank++;
						rankList.put(rank,consensusIFileName);
					}
					else 
						if((!rankList.containsValue(consensusIFileName)) && (!rankList.containsValue(recalibrateJFileName)))
						{
							continue;
						}
			}
		}
		System.out.println("\nPrinting RankList");
		for(Map.Entry<Integer, String> item : rankList.entrySet())
			System.out.println(item);
		System.out.print("\n");
		return rankList;
	}


	/**
	 * Generates a consensus spectrum when two spectra fileNames are provided 
	 * 
	 * @param recalibrateJFileName
	 * @param consensusIFileName
	 * @return SimpleSpectrum
	 */
	public SimpleSpectrum generateConsensus(String recalibrateJFileName,
			String consensusIFileName) {
		RawDataStorage dsrConsensus;
		RawDataStorage dsrRecalibrate;

		dsrConsensus = id.extractMzIntensity(consensusIFileName);
		float[] mzFloatListConsensus = extractFloatMzList(dsrConsensus);
		double[] mzDoubleListConsensus = new double[mzFloatListConsensus.length];
		for(int i = 0; i< mzFloatListConsensus.length;i++)
			mzDoubleListConsensus[i] = mzFloatListConsensus[i];

		float[] intensityFloatListConsensus = extractFloatIntList(dsrConsensus);
		double[] intensityDoubleListConsensus = new double[intensityFloatListConsensus.length];
		for(int i = 0; i< intensityFloatListConsensus.length;i++)
			intensityDoubleListConsensus[i] = intensityFloatListConsensus[i];


		SimpleSpectrum ssConsensus = new SimpleSpectrum(mzDoubleListConsensus,intensityDoubleListConsensus);	

		dsrRecalibrate = id.extractMzIntensity(recalibrateJFileName);
		float[] mzFloatListRecalibrate = extractFloatMzList(dsrRecalibrate);
		double[] mzDoubleListRecalibrate = new double[mzFloatListRecalibrate.length];
		for(int i = 0; i< mzFloatListRecalibrate.length;i++)
			mzDoubleListRecalibrate[i] = mzFloatListRecalibrate[i];

		float[] intensityFloatListRecalibrate = extractFloatIntList(dsrRecalibrate);
		double[] intensityDoubleListRecalibrate = new double[intensityFloatListRecalibrate.length];
		for(int i = 0; i< intensityFloatListRecalibrate.length;i++)
			intensityDoubleListRecalibrate[i] = intensityFloatListRecalibrate[i];

		SimpleSpectrum ssRecalibrate = new SimpleSpectrum(mzDoubleListRecalibrate,intensityDoubleListRecalibrate);

		MutableSpectrum<Peak> newSMutable = performRecalibration(ssConsensus,
				mzDoubleListRecalibrate, intensityDoubleListRecalibrate,
				ssRecalibrate);	

		SimpleSpectrum newConsensus;
		newConsensus = mergeSpectrum(ssConsensus, newSMutable);

		System.out.println("Merged Spectrum generated for " + consensusIFileName + " and " + recalibrateJFileName);
		return newConsensus;
	}


	/**
	 * Generates a consensus spectrum when two spectra are provided 
	 * 
	 * @param newConsensus
	 * @param recalibrateJFileName
	 * @return 
	 */
	public SimpleSpectrum generateConsensus(SimpleSpectrum newConsensus, String recalibrateJFileName) {

		RawDataStorage dsrRecalibrate;

		dsrRecalibrate = id.extractMzIntensity(recalibrateJFileName);
		float[] mzFloatListRecalibrate = extractFloatMzList(dsrRecalibrate);
		double[] mzDoubleListRecalibrate = new double[mzFloatListRecalibrate.length];
		for(int i = 0; i< mzFloatListRecalibrate.length;i++)
			mzDoubleListRecalibrate[i] = mzFloatListRecalibrate[i];

		float[] intensityFloatListRecalibrate = extractFloatIntList(dsrRecalibrate);
		double[] intensityDoubleListRecalibrate = new double[intensityFloatListRecalibrate.length];
		for(int i = 0; i< intensityFloatListRecalibrate.length;i++)
			intensityDoubleListRecalibrate[i] = intensityFloatListRecalibrate[i];

		SimpleSpectrum ssRecalibrate = new SimpleSpectrum(mzDoubleListRecalibrate,intensityDoubleListRecalibrate);

		MutableSpectrum<Peak> newSMutable = performRecalibration(newConsensus,
				mzDoubleListRecalibrate, intensityDoubleListRecalibrate,
				ssRecalibrate);

		newConsensus = mergeSpectrum(newConsensus, newSMutable);

		System.out.println("Merged Spectrum generated for newConsensus and " + recalibrateJFileName);
		return newConsensus;
	}


	/**
	 * Performs recalibration steps when 2 spectra are provided
	 * 
	 * @param ssConsensus
	 * @param mzDoubleListRecalibrate
	 * @param intensityDoubleListRecalibrate
	 * @param ssRecalibrate
	 * 
	 * @return MutableSpectrum
	 */
	public MutableSpectrum<Peak> performRecalibration(
			SimpleSpectrum ssConsensus, double[] mzDoubleListRecalibrate,
			double[] intensityDoubleListRecalibrate,
			SimpleSpectrum ssRecalibrate) {
		double[][] subsetArray = MzRecalibration.maxLinePairStabbing(ssRecalibrate, ssConsensus, 0.5);
		double[] subsetX = new double[subsetArray[0].length];
		double[] subsetY = new double[subsetArray[0].length];

		for(int i=0;i<subsetArray.length;i++)
		{
			if(i==0)
			{
				for(int j=0;j<subsetArray[i].length;j++)
				{
					subsetX[j] = subsetArray[i][j];
				}
			}
			else
				if(i==1)
				{
					for(int j=0;j<subsetArray[i].length;j++)
					{
						subsetY[j] = subsetArray[i][j];
					}
				}
		}

		UnivariateFunction uf = MzRecalibration.getLinearRecalibration(subsetX,subsetY);
		SimpleMutableSpectrum sMutable = new SimpleMutableSpectrum(mzDoubleListRecalibrate,intensityDoubleListRecalibrate);

		MutableSpectrum<Peak> newSMutable =  MzRecalibration.recalibrate(sMutable, uf);
		return newSMutable;
	}


	/**
	 * Merges two spectra (as of here: old consensus and new consensus) 
	 * 
	 * @param SimpleSpectrum oldConsensus, MutableSpectrum<Peak> mutatedSpectrum
	 * @return SimpleSpectrum
	 */		

	public SimpleSpectrum mergeSpectrum(SimpleSpectrum oldConsensus, MutableSpectrum<Peak> mutatedSpectrum)
	{
		ArrayList<Float> mzListNewConsensus = new ArrayList<Float>(); 
		ArrayList<Float> intListNewConsensus = new ArrayList<Float>();
		ArrayList<Float> mzTempStorageList = new ArrayList<Float>();
		ArrayList<Float> intTempStorageList = new ArrayList<Float>();

		TraverseSparseSet tss = new TraverseSparseSet();

		float[][] mzIntOldConsensus;
		float[][] mzIntmutatedSpectrum;

		mzIntOldConsensus = tss.getMzIntensityFromObject(oldConsensus);
		mzIntmutatedSpectrum = tss.getMzIntensityFromObject(mutatedSpectrum);


		float[] mzArrayOld = new float[oldConsensus.size()];	
		float[] intensityArrayOld = new float[oldConsensus.size()];

		for(int i=0;i<mzIntOldConsensus.length;i++)
		{
			for(int j=0;j<mzIntOldConsensus[i].length;j++)
			{
				if (j==0)
				{
					mzArrayOld[i] = mzIntOldConsensus[i][j];
				}
				else
					if(j==1)
					{
						intensityArrayOld[i] = mzIntOldConsensus[i][j];
					}
			}
		}

		float[] mzArrayMutated = new float[mutatedSpectrum.size()];	
		float[] intensityArrayMutated = new float[mutatedSpectrum.size()];


		for(int i=0;i<mzIntmutatedSpectrum.length;i++)
		{
			for(int j=0;j<mzIntmutatedSpectrum[i].length;j++)
			{
				if (j==0)
				{
					mzArrayMutated[i] = mzIntmutatedSpectrum[i][j];
				}
				else
					if(j==1)
					{
						intensityArrayMutated[i] = mzIntmutatedSpectrum[i][j];
					}
			}
		}
		

		//Note: At this point we have four arrays:
		// mzArrayOld, intensityArrayOld, mzArrayMutated, intensityArrayMutated
		double epsilon = 0.5;
		for(int i = 0; i < mzArrayOld.length;i++)
		{
			for(int j = 0; j<mzArrayMutated.length;j++)
			{
				//				if(mzArrayOld[i] == mzArrayMutated[j] || (Math.abs(mzArrayOld[i] - mzArrayMutated[j]) <= epsilon))			
				//				{
				//					mzTempStorageList.add(mzArrayOld[i]);
				//					mzTempStorageList.add(mzArrayMutated[j]);
				//					intTempStorageList.add(intensityArrayOld[i]);
				//					intTempStorageList.add(intensityArrayMutated[j]);
				//
				//					mzListNewConsensus.add(((intensityArrayOld[i] * mzArrayOld[i])+ (intensityArrayMutated[j] * mzArrayMutated[j]))/(intensityArrayOld[i] +intensityArrayMutated[j]));
				//					intListNewConsensus.add(intensityArrayOld[i] +intensityArrayMutated[j]);
				//				}
				//				else if(((mzArrayOld[i] < mzArrayMutated[j])) && (Math.abs(mzArrayOld[i] - mzArrayMutated[j])> epsilon))				
				//				{
				//					mzListNewConsensus.add(mzArrayOld[i]);
				//					intListNewConsensus.add(intensityArrayOld[i]);
				//
				//					mzListNewConsensus.add(mzArrayMutated[j]);
				//					intListNewConsensus.add(intensityArrayMutated[j]);
				//					break;
				//				}
				//				else
				//				{
				//					mzListNewConsensus.add(mzArrayOld[i]);
				//					intListNewConsensus.add(intensityArrayOld[i]);
				//
				//					mzListNewConsensus.add(mzArrayMutated[j]);
				//					intListNewConsensus.add(intensityArrayMutated[j]);
				//					continue;
				//				}
				
				//Condition 1
				if(mzArrayOld[i] < mzArrayMutated[j] || (Math.abs(mzArrayOld[i] - mzArrayMutated[j]) > epsilon))
				{
					mzListNewConsensus.add(mzArrayOld[i]);
					intListNewConsensus.add(intensityArrayOld[i]);	

					mzListNewConsensus.add(mzArrayMutated[j]);
					intListNewConsensus.add(intensityArrayMutated[j]);
					continue;
				}
				else 
					//Condition 2
					if(Math.abs(mzArrayOld[i] - mzArrayMutated[j]) <= epsilon)
					{
						mzTempStorageList.add(mzArrayOld[i]);
						mzTempStorageList.add(mzArrayMutated[j]);

						intTempStorageList.add(intensityArrayOld[i]);
						intTempStorageList.add(intensityArrayMutated[j]);
						
					//	System.out.println("mzArrayOld[i] " + mzArrayOld[i] + " mzArrayMutated[j] " + mzArrayMutated[j]);

						mzListNewConsensus.add(((intensityArrayOld[i] * mzArrayOld[i])+ (intensityArrayMutated[j] * mzArrayMutated[j]))/(intensityArrayOld[i] +intensityArrayMutated[j]));
						intListNewConsensus.add(intensityArrayOld[i] +intensityArrayMutated[j]);
						break;
					}
					else
						//Condition 3
						if(mzArrayOld[i] > mzArrayMutated[j] || (Math.abs(mzArrayOld[i] - mzArrayMutated[j]) > epsilon))
						{
							mzListNewConsensus.add(mzArrayMutated[j]);
							intListNewConsensus.add(intensityArrayMutated[j]);
							break;
						}		
			}
		}

		//		for(int k=0; k<(mzListNewConsensus.size()-1); k++)
		//		{
		//			// Calculate the difference between matchList value k and k+1
		//
		//			float sub1 = mzListNewConsensus.get(k)-mzListNewConsensus.get(k+1);
		//			float sub2 = intListNewConsensus.get(k)-intListNewConsensus.get(k+1);
		//
		//			// If difference is not equal to zero - means the two values are not equal
		//			if((sub1 != 0) && (sub2 != 0) )
		//			{
		//				continue; // iterate to the next value of k
		//
		//			}
		//			else
		//				// If the two values are equal
		//			{	
		//				mzListNewConsensus.remove(k+1); 
		//				intListNewConsensus.remove(k+1);
		//			}
		//		}		
		
		ArrayList<String> stringStorage = new ArrayList<String>(); 
		for (int l = 0;l < mzListNewConsensus.size();l++)
		{
			stringStorage.add(mzListNewConsensus.get(l) + "\t" + intListNewConsensus.get(l));
		}


		mzListNewConsensus.removeAll(mzListNewConsensus);
		intListNewConsensus.removeAll(intListNewConsensus);

		stringStorage = removeDuplicate(stringStorage);
		Collections.sort(stringStorage);

		for(String item : stringStorage)
		{
			String delimiter = "\\t";
			String var[] = item.split(delimiter);		
			String var1 = var[0];
			String var2 = var[1];
			mzListNewConsensus.add(Float.parseFloat(var1));
			intListNewConsensus.add(Float.parseFloat(var2));		
		}
		

		counter=mzListNewConsensus.size();

//		for(float item: mzTempStorageList)
//			System.out.println(item);

		for(int l=0; l<mzTempStorageList.size(); l++)
		{
			for(int k=0; k<counter; k++)
			{
				if(mzListNewConsensus.get(k).equals(mzTempStorageList.get(l)))
				{
					mzListNewConsensus.remove(k);			
					intListNewConsensus.remove(k);
					counter = mzListNewConsensus.size();
				}
			}
		}

		double[] mzArrayNewConsensus = new double[mzListNewConsensus.toArray().length];	
		double[] intensityArrayNewConsensus = new double[mzListNewConsensus.toArray().length];

		for(int i =0; i< mzListNewConsensus.size(); i++)
		{
			mzArrayNewConsensus[i] = mzListNewConsensus.get(i);
			intensityArrayNewConsensus[i] = intListNewConsensus.get(i);
		}

		mzTempStorageList = null;
		intTempStorageList = null;

		SimpleSpectrum newConsensus = new SimpleSpectrum(mzArrayNewConsensus, intensityArrayNewConsensus);
//		System.out.println("Size of new consensus :" + newConsensus.size() + " Printing it now");
//		for(float item : mzListNewConsensus)
//			System.out.println(item);

		return newConsensus;
	}


	/**
	 * Removes duplicate values in a String ArrayList
	 * 
	 * @param ArrayList<String> arlList
	 * @return ArrayList<String>
	 */	
	public ArrayList<String> removeDuplicate(ArrayList<String> arlList)
	{
		HashSet<String> h = new HashSet<String>(arlList);
		arlList.clear();
		arlList.addAll(h);
		return arlList;
	}


	/**
	 * Extracts the mz and intensity values from a MutableSpectrum<Peak> object  
	 * 
	 * @param MutableSpectrum<Peak> 
	 * @return double[][] mzIntensityArray
	 */	
	public float[][] getMzIntensityFromObject(MutableSpectrum<Peak> ms)
	{
		float[][] mzIntensityArray = new float[ms.size()][2];
		for(int i=0; i< ms.size();i++)
		{
			mzIntensityArray[i][0] = (float) ms.getMzAt(i);
			mzIntensityArray[i][1] = (float) ms.getIntensityAt(i);
		}
		return mzIntensityArray;
	} 


	/**
	 * Extracts the mz and intensity values from a SimpleSpectrum object 
	 * 
	 * @param SimpleSpectrum 
	 * @return double[][] mzIntensityArray
	 */	
	public float[][] getMzIntensityFromObject(SimpleSpectrum ss)
	{
		float[][] mzIntensityArray = new float[ss.size()][2];
		for(int i=0; i< ss.size();i++)
		{
			mzIntensityArray[i][0] = (float) ss.getMzAt(i);
			mzIntensityArray[i][1] = (float) ss.getIntensityAt(i);
		}
		return mzIntensityArray;
	}


	/**
	 * Extracts the double intensity array from the RawDataStorage object 
	 * 
	 * @param RawDataStorage dsr
	 * @return double[] intensityDoubleList
	 */	

	public float[] extractFloatIntList(RawDataStorage dsr) {
		float[] intensityFloatList = new float[dsr.mzArray.length];

		int j=0;
		for(float arrayItem : dsr.intensityArray){
			intensityFloatList[j] = arrayItem;
			j++;
		}
		return intensityFloatList;
	}


	/**
	 * Extracts the double mz array from the RawDataStorage object 
	 * 
	 * @param RawDataStorage dsr
	 * @return double[] mzDoubleList
	 */	

	public float[] extractFloatMzList(RawDataStorage dsr) {
		float[] mzFloatList = new float[dsr.mzArray.length];

		int d=0;
		for(float arrayItem : dsr.mzArray){
			mzFloatList[d] = arrayItem;
			d++;
		}	
		return mzFloatList;
	}
}
