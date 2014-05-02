/**
 * 
 */

package de.uni_jena.bio.informatik.input;

import java.io.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import de.uni_jena.bio.informatik.functions.CalculateDistance;
import de.uni_jena.bio.informatik.functions.MatrixFunctions;
import de.uni_jena.bio.informatik.functions.SparseMatrix;
import de.uni_jena.bio.informatik.functions.TraverseSparseSet;
import de.unijena.bioinf.ChemistryBase.ms.MutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;


/**
 * InputData.
 * <h3>Usage</h3>
 * <ol>
 * <li>Accepts user input in the form of directory path {@code userInput}.
 *  Sorts the files in the directory </li>
 * <li>Creates an InputData object {@code InputData} and calls the userInput {@code userInput} 
 * to ask for user input</li>
 * <li>Extracts mz and intensity values as array {@code extractMzIntensity} and checks for empty spectra</li>
 *  <li>Extracts mz and intensity values as array {@code extractMzIntensity}</li>
 *  <li> Generates match list using 2 spectra as input </li>
 *  <li>Calculates distance between two mz lists {@code specDistance}</li>
 *  <li>Transposes the distance matrix {@code transposeFloat}</li>
 *  <li>Generates a sparse matrix of this distance matrix {@code generateSparseMatrix}</li>
 *  <li>Sparse matrix is traversed to perform recalibration
 *  starting from the two mz lists having the minimum distance {@code traverseSortedSet}</li>
 *
 * </ol>
 * 
 * @author Purva Kulkarni
 *
 */

public class InputData {

	static String dirPath;

	public void userInput() throws IOException
	{

		//For user input
//		Scanner scanner = new Scanner(System.in );
//		System.out.println("Enter the file path: ");
//		dirPath = scanner.nextLine(); // Takes the directory path as the user input
		dirPath = "/home/purva/Lab_proceedings/Work/temp";

		File folder = new File(dirPath);
		//checks whether the specified path has a folder
		if(folder.isDirectory())
		{
			File[] fileList = folder.listFiles();
			
			SortFile sf = new SortFile();

			fileList = sf.sortByNumber(fileList);

			ArrayList<String> fileNameList = new ArrayList<String>();
			ArrayList<String> coordinateList = new ArrayList<String>();
			ArrayList<Object> spectrumArrayList = new ArrayList<Object>();
			ArrayList<String> fileDataI = new ArrayList<String>();

			RawDataStorage dsr;
			CalculateDistance cd1 = new CalculateDistance();
			int fileCount = 0;
			for(int i=0;i<(fileList.length);i++)
			{	
				File spectrumFile = fileList[i];

				String fileNameI = new String();
				fileNameI = folder + "/" + spectrumFile.getName();
				// System.out.println("Filename"+ i+":"+spectrumFile.getName());

				ReadFile rf = new ReadFile();

				// call the fileInput method and passes fileName to it
				fileDataI = rf.fileInput(fileNameI); //returned ArrayList stored in fileData


				// Call the listTo Array Method
				String[] fileDataArrayI = listToArray(fileDataI);

				float[] mzArray = new float[fileDataArrayI.length];
				float[] intensityArray = new float[fileDataArrayI.length];
				float[] mzList = new float[fileDataArrayI.length];

				//Note: coordinate and count removed as of now	
				dsr = extractMzIntensity(fileDataArrayI, mzArray, intensityArray, spectrumFile.getName());
				mzList = dsr.mzArray;

				if((cd1.emptySpectra(mzList))== false)
				{
					fileCount++;
					continue;
				}
				else
				{
					// Add the mzLists and its corresponding coordinate values to the respective Array Lists
					spectrumArrayList.add(mzList);
					coordinateList.add(dsr.coordinates);
					fileNameList.add(dsr.fileName);
				}

			}
			

			if(fileCount == fileList.length)
				System.out.println("All "+ fileCount + " files contain empty spectra. The program will now terminate!");


			float[][] distanceMatrix = new float[spectrumArrayList.size()][spectrumArrayList.size()];
			String[][] coordinateMatrix = new String[spectrumArrayList.size()][spectrumArrayList.size()];
			MatchListDataStorage mdsObject;


			for(int i=0; i < spectrumArrayList.size();i++)
			{
				Object mzObject = spectrumArrayList.get(i);
				float[] consensusI = (float[])mzObject;

				String coordinateOfI = coordinateList.get(i);
				String FileNameOfI = fileNameList.get(i);

				for(int j=i+1; j < spectrumArrayList.size();j++)
				{
					mzObject = spectrumArrayList.get(j);
					float[] referenceJ = (float[])mzObject;

					String coordinateOfJ = coordinateList.get(j);
					String FileNameOfJ = fileNameList.get(j);

					mdsObject= cd1.specDistance(consensusI,referenceJ, coordinateOfI, coordinateOfJ, FileNameOfI,FileNameOfJ);
					distanceMatrix[i][j] = mdsObject.distanceValue;
					coordinateMatrix[i][j] = mdsObject.twoFileCoordinates;				
				}
			}


			MatrixFunctions tm = new MatrixFunctions();

			// Matrix obtained after transposing
			distanceMatrix = tm.transposeFloat(distanceMatrix, distanceMatrix.length);


			// Coordinate Matrix after transposing
			coordinateMatrix = tm.transposeString(coordinateMatrix, coordinateMatrix.length);

			SparseMatrix sm = new SparseMatrix();
			Set<Map.Entry<String, Float>> sortedsparse = sm.generateSparseMatrix(distanceMatrix,coordinateMatrix);
	//		System.out.println("Printing Sorted Sparse Set");

			PrintWriter out = new PrintWriter(new FileWriter("output.xls"));
		//	PrintWriter out = new PrintWriter(new FileWriter("/home/purvakulkarni/Desktop/output.xls"));

			for(Map.Entry<String, Float> element : sortedsparse)
				out.println(element);

			out.close();

			TraverseSparseSet tss = new TraverseSparseSet();
			SimpleSpectrum finalConsensus = tss.traverseSortedSet(sortedsparse);

			System.out.println("\nFinal Consensus Generated !\n===============================\nNow Recalibrating\n");


			for(String item: fileNameList)
			{
				RawDataStorage dsrObject = extractMzIntensity(item);
				float[] mzArrayRecalibrate = dsrObject.mzArray;
				double[] mzDoubleArrayRecalibrate = new double[mzArrayRecalibrate.length];
				for(int j = 0; j< mzDoubleArrayRecalibrate.length;j++)
					mzDoubleArrayRecalibrate[j] = mzArrayRecalibrate[j];

				float[] intArrayRecalibrate = dsrObject.intensityArray;
				double[] intDoubleArrayRecalibrate = new double[intArrayRecalibrate.length];
				for(int j = 0; j< intDoubleArrayRecalibrate.length;j++)
					intDoubleArrayRecalibrate[j] = intArrayRecalibrate[j];

				SimpleSpectrum ssToRecalibrate = new SimpleSpectrum(mzDoubleArrayRecalibrate,intDoubleArrayRecalibrate);

				TraverseSparseSet tssObject = new TraverseSparseSet();

				MutableSpectrum<Peak> newRecalibratedSpectrum = tssObject.performRecalibration(finalConsensus,
						mzDoubleArrayRecalibrate, intDoubleArrayRecalibrate,
						ssToRecalibrate);	

				double[] newRecalibratedSpectraMz = new double[newRecalibratedSpectrum.size()];
				double[] newRecalibratedSpectraInt = new double[newRecalibratedSpectrum.size()];

				float[] newRecalibratedSpectraMzFloat = new float[newRecalibratedSpectrum.size()];
				float[] newRecalibratedSpectraIntFloat = new float[newRecalibratedSpectrum.size()];


				for(int k = 0; k<newRecalibratedSpectrum.size();k++)
				{
					newRecalibratedSpectraMz[k] = newRecalibratedSpectrum.getMzAt(k);
					newRecalibratedSpectraMzFloat[k] = (float) newRecalibratedSpectraMz[k];

					newRecalibratedSpectraInt[k] = newRecalibratedSpectrum.getIntensityAt(k);
					newRecalibratedSpectraIntFloat[k] = (float) newRecalibratedSpectraInt[k];
				}

				File dir = new File("RecalibratedSpectra/");
				dir.mkdir();
				System.out.println("New folder created");
				
				PrintWriter out1 = new PrintWriter(new FileWriter( dir +"/"+ item));

				for(int l = 0; l < newRecalibratedSpectraMzFloat.length;l++)
					out1.println(newRecalibratedSpectraMzFloat[l]+ "\t\t" + newRecalibratedSpectraIntFloat[l]);		
				out1.close();
			}

			System.out.println("\nAll new recalibrated files generated");


		}
	}


	/**
	 * Converts the Array List to array and returns a string array (mz and intensity)
	 * 
	 * @param ArrayList<String> fileData
	 * @return fileDataArray 
	 */	

	private String[] listToArray(ArrayList<String> fileData) {
		//Creates a new Array fileDataArray
		String[] fileDataArray = new String[fileData.size()];
		fileDataArray = fileData.toArray(fileDataArray); //Converts ArrayList to Array
		return fileDataArray;	
	}   


	/**
	 * Extracts m/z and intensity list - method overriding 1.
	 * Returns object of class RawDataStorage {@link RawDataStorage} containing the following fields
	 * coordinateString, newMzArray, newIntensityArray, FileName.
	 * 
	 * @param String[] fileDataArray
	 * @param float[] mzArray
	 * @param float[] intensityArray
	 * @param String FileName
	 * @return RawDataStorage 
	 */	

	private RawDataStorage extractMzIntensity(String[] fileDataArray, float[] mzArray,
			float[] intensityArray, String FileName) {

		float[] newMzArray = new float[mzArray.length-3];
		float[] newIntensityArray = new float[intensityArray.length-3];
		String coordinateString = null;

		//	System.out.println("File:" + FileName );
		for(int temp=0; temp<fileDataArray.length;temp++)
		{
			// If the line count is between 1 and 3, split the variables as x and y
			if((temp > 0) && (temp < 2)) {	
				coordinateString = fileDataArray[temp];		
			}	
			// If the line count is greater than 3, print the lines
			else
				if(temp > 2) 
				{
					String var[] = fileDataArray[temp].split("		");		
					String var3 = var[0];
					String var4 = var[1];
					mzArray[temp] = Float.parseFloat(var3);
					intensityArray[temp] = Float.parseFloat(var4);
				}
		}
		int tempVal=3;
		for(int l=0; l < newMzArray.length;l++)
		{
			newMzArray[l]= mzArray[tempVal];
			newIntensityArray[l]= intensityArray[tempVal];
			tempVal++;
		}

		RawDataStorage ds = new RawDataStorage(coordinateString, newMzArray, newIntensityArray, FileName);
		return ds;
	}


	/**
	 * Extracts m/z and intensity list when only a single file name is provided - method overriding 2.
	 * Returns object of class RawDataStorage {@link RawDataStorage} containing the following fields
	 * fileDataArrayI, mzArray, intensityArray, fileName.
	 * 
	 * @param fileName
	 * @return RawDataStorage 
	 */	

	public RawDataStorage extractMzIntensity(String fileName) {

		String fullFilePath = dirPath + "/" + fileName;
		ArrayList<String> fileDataI = new ArrayList<String>();
		RawDataStorage dsr;

		ReadFile rfI = new ReadFile();

		// call the fileInput method and passes fileName to it
		fileDataI = rfI.fileInput(fullFilePath); //returned ArrayList stored in fileData

		// Call the listTo Array Method
		String[] fileDataArrayI = listToArray(fileDataI);

		float[] mzArray = new float[fileDataArrayI.length];
		float[] intensityArray = new float[fileDataArrayI.length];

		dsr = extractMzIntensity(fileDataArrayI, mzArray, intensityArray, fileName);
		return dsr;
	}	
}
