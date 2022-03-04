/*E.Camerini 06/2020 version 4.0
 * Addapted by M. Svoren 01/2022 to version 4.1
 * 
This macro yields a table containing several 
nuclear shape descriptors and the Nuclear irregulatiry index.
run this macro once the raw image is open.*/ 

//BEFORE STARTING CHANGE THE NUMBER IN "C*-AVG_" TO MATCH THE CHANNEL TO THE SIGNAL

			Nucleus= 4; //use nuclear stain signal 
			
			
//ask for file path directory
			
			imageTitle=getTitle();//returns a string with the image title
			barename=replace(imageTitle, ".czi",""); //used to save results
			getPixelSize(unit, pixelWidth, pixelHeight);//needed for distance calculations

//Z stack and separate channels
			run("Z Project...", "projection=[Average Intensity]");
			run("Duplicate...", "duplicate channels=[Nucleus]");
			rename("DAPI");

// Threshold the cells shapes. Li was manually chosen before applied to all samples
			run("Median...", "radius=1 stack");
			setAutoThreshold("Li dark");
			call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			run("Duplicate...", " ");
			rename("DAPI threshold");
			//run("Watershed"); for this set of samples the watershed was removed from the analysis, 
			// as a consequence the unsegmented cells need to be manually deselected
			run("Set Measurements...", "area mean min center bounding shape feret's redirect=None decimal=3");
			run("Analyze Particles...", "size=50-Infinity display exclude clear add");

// here the macro stops to give you time to delete the ROIs of unsegmented cells manually, press ok when done
			demo = "WaitForUserDemo";
			msg = "If necessary, deselect unsegmented groups of cells, then click \"OK\".";
			waitForUser(demo, msg);
			outPath = getDirectory("Choose a Directory to save Results");
			roiManager("save", outPath + barename + "_RoiSet.zip");
			setBatchMode(true);
			run("Clear Results");
			selectWindow("DAPI") 
			run("Set Measurements...", "area mean min center bounding shape feret's redirect=None decimal=3");
			roiManager("Deselect");
			roiManager("Measure");
			IJ.renameResults("Original set");

// Going through all the nuclei and enlarge each ROI by 1 pixel, then adding this enlarged ROI to the manager

	selectWindow("DAPI")
	nROIs=roiManager("count");
	print(nROIs);

	for(i=0; i<nROIs; i++){
		roiManager("Select", i);
		run("Enlarge...", "enlarge=1 pixel");
		roiManager("Add");
	}

//selecting the enlarged ROIS 

	nROIs2=roiManager("count");
	for(i=nROIs; i<nROIs2; i++){
		roiManager("Select", i);
		run("Set Measurements...", "area mean min center bounding shape feret's redirect=None decimal=3");
		roiManager("Measure");
	}
	
	selectWindow("Results");
	Table.rename("Results", "Enlarged Results");


//In order to calculate the smallest and largest radius in the cell a 1 pixel wide ring around each nucleus is selected by subtracting
//the original selection to the 1-pixel enlarged one. All the pixels in this selection are put in arrays with Roi.getContained Points.
//Then each is subtracted to the x and y value of the centroid coordinates of the nucleus and the smallest and largest are selected.

	for(i=nROIs; i<nROIs2; i++){
		roiManager("Select", i);
		run("Create Mask");
		imageCalculator("Subtract", "Mask", "DAPI threshold"); 
		run("Create Selection"); //Creating and adding the ring selection to the ROI manager
		roiManager("Add");
		
		Roi.getContainedPoints(xpoints, ypoints);
		maxDist=0;
		minDist=10000000000000000000000; //indicative number to ensure the radius is smaller than this
	
		selectWindow("Original set");
		Xcentroid= parseFloat(Table.get("XM",i-nROIs));
		Ycentroid= parseFloat(Table.get("YM",i-nROIs));
		
			for (j = 0; j < xpoints.length; j++) {
				xDiff= Xcentroid - (xpoints[j]* pixelWidth);
				yDiff= Ycentroid - (ypoints[j]* pixelHeight);
				dist= sqrt((xDiff*xDiff)+(yDiff*yDiff));
				if(maxDist< dist){
					maxDist=dist;
				}
					if(minDist> dist){
						minDist=dist;
					}
			}

			//write the smallest and largest radii of each nuclei into a result table
			setResult("RefNumber", i-nROIs, i-nROIs+1);
			setResult("MaxDist", i-nROIs, maxDist);
			setResult("MinDist", i-nROIs, minDist);
			updateResults();
			
			//to avoid a thousand windows opening
			if (i<nROIs2-1) {
				selectWindow("Mask");
				close(); //close all but the last
			}
	}

	selectWindow("Results");
	Table.rename("Results", "Cell Radii");

//Combining the results needed in a final table. Each column starts from an array that loops through 0->nROIs

	area=newArray(nROIs);
	aspectRatio=newArray(nROIs);
	radiusMax=newArray(nROIs);
	radiusMin=newArray(nROIs);
	radiiRatio=newArray(nROIs);
	circularity=newArray(nROIs);
	boxWidth=newArray(nROIs);
	boxHeight=newArray(nROIs);
	roundness=newArray(nROIs);
	solidity=newArray(nROIs);

	for(k=0; k<nROIs; k++){
		selectWindow("Original set");
		area[k]= parseFloat(Table.get("Area",k));
		aspectRatio[k]= parseFloat(Table.get("AR",k));
		circularity[k]= parseFloat(Table.get("Circ.",k));
		boxWidth[k]= parseFloat(Table.get("Width",k));
		boxHeight[k]= parseFloat(Table.get("Height",k));
		roundness[k]= parseFloat(Table.get("Round",k));
		solidity[k]= parseFloat(Table.get("Solidity",k));
		selectWindow("Cell Radii");
		radiusMax[k]= parseFloat(Table.get("MaxDist",k));
		radiusMin[k]= parseFloat(Table.get("MinDist",k));
		radiiRatio[k]=radiusMax[k]/radiusMin[k];
	}

Table.create("Final results");
	for(k=0; k< nROIs; k++){
		Table.set("Refnumb", k, k+1);
		Table.set("Area", k, area[k]);
		Table.set("Aspect ratio", k, aspectRatio[k]);
		Table.set("Radii ratio", k, radiiRatio[k]);
		Table.set("Circularity^-1", k, 1/circularity[k]);
		Table.set("Area over box", k, area[k]/(boxWidth[k]*boxHeight[k]));
		Table.set("NII",k,aspectRatio[k]+radiiRatio[k]+(1/circularity[k])-(area[k]/(boxWidth[k]*boxHeight[k])));
		Table.set("Roundness", k, roundness[k]);
		Table.set("Solidity", k, solidity[k]);
	}
	
	setBatchMode(false);
	//Done!
	selectWindow("Final results");
	saveAs("Results", outPath + barename +"_CompleteResults.txt");
		
		selectWindow(barename +"_CompleteResults.txt");	
		run("Close");
		     	 
         selectWindow("Cell Radii"); 
         run("Close");
		    	
         selectWindow("Enlarged Results"); 
         run("Close");
    	    	 
         selectWindow("Original set"); 
         run("Close");
				
         selectWindow("ROI Manager");
         run("Close");
    	   		
         selectWindow("Log");
         run("Close");

         run("Close All");
         
        
		 
