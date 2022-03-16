/*
This macro can be used to analyse nucleus shape change over time in a sequential image stack.
The input is a thresholded mask of the original image stack and a labeled, point-wise tracking stack (as created by, for example, the TrackMate plugin).
The tracking stack of the image sequence is used to identify and follow each nuclear segment in the subsequent time frames. 
*/

// Switch USE_DEFAULT_NAMES to true in order to use standard names for the segmented mask and the tracked spots stacks.
// This no longer requires the user to select the mask and tracking image by hand.
// You can change the MASK_TITLE and SPOTS_TITLE to match your own standard image names.
USE_DEFAULT_NAMES = false;
MASK_TITLE = "MaskStack.tif";
SPOTS_TITLE = "TrackedSpots.tif";

// Settings for the mask creation
CREATE_MASK = true; // Set to true if the chosen mask image still needs to be thresholded. The next few settings will determine how. If set to false, the mask image is assumed to already be a binary image and the thresholding will be skipped.
PREPRO_FILTER = "Median..."; // This is the preprocessing filter command. Normally "Median..." works best, but here you can change it if need be.
FILTER_RADIUS = 5; // The filter radius (in pixels)
THRESHOLD_METHOD = "Otsu"; // The name of the threshold used. These names can be found in the dialog of the Image > Adjust > Threshold... command.

STORE_TRACKED_ROIS = true; // Set to true in order to save the outlines of the tracked segments per track, with roi names reflecting the track id and frame number of the segment. Set to false to not store these ROIs.

// Preset constants for the resulting table and its contents. Do not adapt these without also changing the macro to suit.
RESULTS_TABLE = "Nuclear shape change analysis";
TAB_HEADER = newArray("Track ID", "Frame Nr", "Area", "Aspect Ratio", "Radii Ratio", "Circularity^-1", "Area over box", "NII", "Delta NII", "Roundness", "Solidity");

//--------------------------------------------------------------------------------

// Get the proper images based on settings or user interactions
if(USE_DEFAULT_NAMES)
{
	selectWindow(MASK_TITLE);
	maskID = getImageID();
	getPixelSize(unit, pixelWidth, pixelHeight);//needed for distance calculations
	selectWindow(SPOTS_TITLE);
	spotsID = getImageID();
}
else
{
	// Get mask image ID
	Dialog.createNonBlocking("Please select the mask image");
	Dialog.addMessage("Please select the mask image.");
	Dialog.show();
	maskID = getImageID();
	title = getTitle();
	title = File.getNameWithoutExtension(title);
	getPixelSize(unit, pixelWidth, pixelHeight);//needed for distance calculations

	// Get mask image ID
	Dialog.createNonBlocking("Please select the tracked spots image");
	Dialog.addMessage("Please select the tracked spots image.");
	Dialog.show();
	spotsID = getImageID();
}

resultsPath = getDir("Please select a results folder");

setBatchMode(true);
nrOfSlices = nSlices;  // Store the number of slices in the movie

// Create the mask image
if(CREATE_MASK)
{
	selectImage(maskID);
	run(PREPRO_FILTER, "radius=" + FILTER_RADIUS + " stack");
	setOption("BlackBackground", true);
    run("Convert to Mask", "method=" + THRESHOLD_METHOD + " background=Dark calculate black");
}

// Now let's use the tracked spots to segment the mask movie
for(i = 0; i < nrOfSlices; i++)
{
	// Duplicate the current slice of both movies to prepare the segmentation for this time frame
	selectImage(maskID);
	setSlice(i + 1);
	run("Duplicate...", "title=Mask");
	selectImage(spotsID);
	setSlice(i + 1);
	run("Duplicate...", "title=Marker");
	// Segment the mask slice with the tracked spots as seeds
	run("Marker-controlled Watershed", "input=Mask marker=Marker mask=Mask");
	// Set min intensity to 0 as the marker watershed gives non-marked (i.e. untracked) segments a -1 value.
	run("Min...", "value=0");

	// Create the segmented movie
	if(i == 0)
	{
		// First slice, start the new stack
		rename("WatershedResult");
	}
	else 
	{
		// New slice added to existing segmented stack
		curTitle = getTitle();
		run("Concatenate...", "  title=WatershedResult image1=WatershedResult image2=" + curTitle + " image3=[-- None --]");
	}

	// Clean up by closing the copies of the current slice
	selectWindow("Mask");
	close();
	selectWindow("Marker");
	close();
}

// Find the highest label ID in the stack
maxLabel = 0;
for(i = 1; i <= nSlices; i++)
{
	setSlice(i);
	getRawStatistics(nPixels, mean, min, maxSliceLabel, std, histogram);
	if(maxSliceLabel > maxLabel)
	{
		maxLabel = maxSliceLabel;
	}
}

resultsTable = Table.create(RESULTS_TABLE);

// For all labels, do all the measurements in turn
for (i = 1; i <= maxLabel; i++) 
{
	// Select the current label from the labeled watershed results and create a mask out of it
	selectWindow("WatershedResult");
	run("Select Label(s)", "label(s)=" + i);
	setThreshold(0.5000, 1000000000000000000000000000000.0000);
	run("Convert to Mask", "method=Default background=Dark black");

	// Measure the created segments in the entire stack. This creates a results table with one line of measurements for each frame the nucleus appeared in.
	run("Set Measurements...", "frame area mean min center bounding shape feret's stack redirect=None decimal=3");
	run("Analyze Particles...", "clear display stack add");
	if(STORE_TRACKED_ROIS && nResults > 0)
	{
		for(roiNr = 0; roiNr < nResults; roiNr++)
		{
			sliceNr = getResult("Slice", roiNr);
			roiManager("select", roiNr);
			roiManager("rename", "track_" + i + "_frame_" + sliceNr);
		}
		run("Select None");
		roiManager("save", resultsPath + title + "_RoiSet_Track_" + i + ".zip");
	}

	// Reset the 'previous Nii' feature as we start with a new nucleus
	prevNii = NaN;
	// Loop through each measured frame to collect data and calculate the missing features
	for(resNr = 0; resNr < nResults; resNr++)
	{
		// Get the slice value from the results to set the segment stack to the right time frame
		sliceNr = getResult("Slice", resNr);
		setSlice(sliceNr);

		// Get the measured features
		area = getResult("Area", resNr);
		aspectRatio = getResult("AR", resNr);
		circularity = getResult("Circ.", resNr);
		boxWidth = getResult("Width", resNr);
		boxHeight = getResult("Height", resNr);
		roundness = getResult("Round", resNr);
		solidity = getResult("Solidity", resNr);
		xCentroid = getResult("XM", resNr);
		yCentroid = getResult("YM", resNr);

		// Create the 1-pixel-wide border ROI for the current selection
		roiManager("reset");  // Do this explicitly for if no ROIs are found, no clear will have taken place
		// Start by adding the original outline to the roi manager
		run("Create Selection");
		roiManager("Add");
		roiManager("select", 0);
		// Enlarge by one pixel and add this new outline as well
		run("Enlarge...", "enlarge=1 pixel");
		roiManager("Add");
		// Create the XOR of the two outlines. This gives you a one pixel wide selection around the original segment.
		roiManager("Select", newArray(0,1));
		roiManager("XOR");
		// Calculate the min and max radius ratio of this 1-pixel-wide ROI
		radiiRatio = getRadiiRatio(xCentroid, yCentroid, pixelWidth, pixelHeight);

		// Calculate the missing values and add all to the final results table
		prevNii = calculateResults(i, sliceNr, area, aspectRatio, radiiRatio, circularity, boxWidth, boxHeight, roundness, solidity, prevNii);
	}

	close(); // Closes label image
}
close("Results");
selectWindow("WatershedResult");
save(resultsPath + title + "_tracked_segmentation.tif");
close();
selectWindow(RESULTS_TABLE);
Table.save(resultsPath + title + "_" + RESULTS_TABLE + ".csv");

setBatchMode(false);
print("Finished nculear shape change analysis of " + title + ".");

// Function to calculate the ratio between the maximum and minimum radius of a single-pixel-wide circular ROI.
// The radius is calulated from the given coordinates of the center of mass of the internal area enclosed by the ROI.
// @param aXcentroid	The x-coordinate of the center of mass
// @param aYcentroid	The y-coordinate of the center of mass
// @param aPixelWidth	The pixel width of the image in order to correctly calculate the distance between the CoM and the Roi pixels
// @param aPixelHeight	The pixel height of the image in order to correctly calculate the distance between the CoM and the Roi pixels
function getRadiiRatio(aXcentroid, aYcentroid, aPixelWidth, aPixelHeight)
{
	// Get all points from the current ROI. These are all the possible end points for a minimum and maximum radius.
	Roi.getContainedPoints(xpoints, ypoints);

	// Set the starting values of the min and max radius so that they will always be superseded by any calculated radius.
	maxDist=0;
	minDist=10000000000000000000000;
	
	
	for (j = 0; j < xpoints.length; j++) 
	{
		xDiff= aXcentroid - (xpoints[j]* aPixelWidth);
		yDiff= aYcentroid - (ypoints[j]* aPixelHeight);
		dist= sqrt((xDiff*xDiff)+(yDiff*yDiff));
		if(maxDist< dist)
		{
			maxDist=dist;
		}
		if(minDist> dist)
		{
			minDist=dist;
		}
	}

	return maxDist/minDist;
}

// Calculate area over box, the NII and the NII delta based on the measured fetures for one nucleus in ne time-frame. Add all values to the results table.
// @param	aTrackID		The id of the tracked nucleus
// @param	aFrameNr		The frame number that these measurements are taken from
// @param	aArea			The area of the nucleus
// @param	aAR				The aspect ratio (longest vs shortest axis of the best-fit ellipse)
// @param	aRadiiRatio		The ratio between the maximum vs the minimum radius of the nucleus
// @param	aCirc			The circularity measure of the nucleus
// @param	aWidth			The width of the bounding box around the nucleus
// @param	aHeight			The height of the bounding box around the nucleus
// @param	aRoundness		The roundness measure of the nucleus
// @param	aSolidity		The solidity measure of the nucleus
// @param	aPrevNii		The value of the NII of this nucleus in the previous time frame
//
// @return	The calculated Nii value
function calculateResults(aTrackID, aFrameNr, aArea, aAR, aRadiiRatio, aCirc, aWidth, aHeight, aRoundness, aSolidity, aPrevNii)
{
	circInv = 1/aCirc;
	areaOverBox = aArea/(aWidth * aHeight);
	nii = aAR + aRadiiRatio + circInv - areaOverBox;
	// Only calculate delta Nii if there is a previous value
	if(!isNaN(aPrevNii))
	{
		deltaNii = abs(aPrevNii - nii);
	}
	else 
	{
		deltaNii = NaN;
	}
	
	selectWindow(RESULTS_TABLE);
	nextRow = Table.size;
	Table.set(TAB_HEADER[0], nextRow, aTrackID);
	Table.set(TAB_HEADER[1], nextRow, aFrameNr);
	Table.set(TAB_HEADER[2], nextRow, aArea);
	Table.set(TAB_HEADER[3], nextRow, aAR);
	Table.set(TAB_HEADER[4], nextRow, aRadiiRatio);
	Table.set(TAB_HEADER[5], nextRow, circInv);
	Table.set(TAB_HEADER[6], nextRow, areaOverBox);
	Table.set(TAB_HEADER[7], nextRow, nii);
	Table.set(TAB_HEADER[8], nextRow, deltaNii);
	Table.set(TAB_HEADER[9], nextRow, aRoundness);
	Table.set(TAB_HEADER[10], nextRow, aSolidity);

	Table.update;

	return nii; 
}
			