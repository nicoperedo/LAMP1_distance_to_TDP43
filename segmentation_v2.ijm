// @ File(label="File directory", style="directory") dir
// @ File(label="ROI directory", style="directory") roi_dir
// @ String (label="Batch process", choices={"Yes", "No"}, style="listBox") batch
// @ Integer (label="File for test", min=0, max=100, value=5) test_file
// @ Integer (label="ROI number for test", min=0, max=10, value=2) test_roi
// @ String (label="File suffix", choices={".nd2", ".tif"}, style="listBox") suffix

//Parameters
nuclei_channel = 1;
spot_channel = 3;
cluster_channel = 2;
//raw_image = "raw";
binsize_sep = 0.5;
binsize_dist = 100;
binsize_dist_scaled = 0.1;

threshold_cluster = 7000;

pixel_size = 0.0693948;
z_step = 0.5;

filter_nuclei = 1000;//filter in microns
filter_nuclei_voxels = filter_nuclei/pixel_size/pixel_size/z_step;

//Column names
column_names = newArray("Values[um]","Counts[px]","Counts_fraction","Cell_id","Aggregate_ratio","Filename");

//File directory
fileList = getFilesList(dir, suffix);
Array.sort(fileList);
//ROI directory
roiList = getFilesList(roi_dir, ".zip");
Array.sort(roiList);

//Create the different folders with results
File.makeDirectory(dir + "/Analysis");
File.makeDirectory(dir + "/Binaries");

if (batch == "Yes") {
	setBatchMode(true);
}

for (files = 0; files < fileList.length; files++) {
//for (files = test_file; files < test_file+1; files++) {
	//File and ROI names
	raw_image = fileList[files];
	drawn_roi = roiList[files];
	name = getBasename(raw_image, suffix);
	
	//Open image
	run("Bio-Formats Importer", "open=[" + dir + File.separator + raw_image + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");

	//Open ROI Manager and load ROI
	
	roiManager("Open", roi_dir + File.separator + drawn_roi);
		
	//Process
	segment_spots_cells(raw_image,spot_channel);
	segment_cluster(raw_image,cluster_channel,threshold_cluster);
	segment_nuclei(raw_image,nuclei_channel,filter_nuclei_voxels);
	process_roi();
	
	
	count = roiManager("count");
	array_count = Array.getSequence(count);
	roiManager("select", array_count);
	roiManager("Delete");
	close("Results");
	run("Collect Garbage");
	run("Collect Garbage");
	close(raw_image);
	close("spots");
	close("cluster");
	close("spots_raw");
	close("totalcell_binary");
	close("nuclei");
}

//close("*");

//Extract a string from another string at the given input smaller string (eg ".")
function getBasename(filename, SubString){
  dotIndex = indexOf(filename, SubString);
  basename = substring(filename, 0, dotIndex);
  return basename;
}

//Return a file list contain in the directory dir filtered by extension.
function getFilesList(dir, fileExtension) {  
  tmplist=getFileList(dir);
  list = newArray(0);
  imageNr=0;
  for (i=0; i<tmplist.length; i++)
  {
    if (endsWith(tmplist[i], fileExtension)==true || endsWith(tmplist[i], fileExtension + ".roi")==true)
    {
      list[imageNr]=tmplist[i];
      imageNr=imageNr+1;
      //print(tmplist[i]);
    }
  }
  Array.sort(list);
  return list;
}

//Apply voxel calibration and units from an image to another. Assumes the image to calibrate has a single channel and no time dimension
function calibrate_image(calibrated_image, image_to_calibrate){
	selectWindow(calibrated_image);
	
	//Get parameters from image
	getDimensions(width, height, channels, slices, frames);
	getVoxelSize(width, height, depth, unit);
	print(unit);
	
	selectWindow(image_to_calibrate);
	
	//Apply calibration
	//Reassign pixel size since the generated filtered image is not calibrated anymore
	Stack.setXUnit(unit);
	run("Properties...", "channels=1 slices=" + slices + " frames=1 pixel_width=" + width + " pixel_height=" + width + " voxel_depth=" + depth);
}

//Segment nuclei
function segment_nuclei(raw_image,nuclei_channel,filter_nuclei_voxels) { 
	//Isolate channel
	selectWindow(raw_image);
	run("Duplicate...", "duplicate channels=" + nuclei_channel + "-" + nuclei_channel); // In the case of Jacqueline 
	rename("nuclei_raw");
	
	//Normalize the entire stack
	run("Enhance Contrast...", "saturated=0 normalize process_all use"); 
	
	//Segment
	run("Gaussian Blur...", "sigma=3 scaled stack");
	setAutoThreshold("Default dark no-reset stack");
	run("Convert to Mask", "background=Dark black");
	
	//fill
	run("Fill Holes", "stack");
	
	//Generate labelled nuclei
	run("Distance Transform Watershed 3D", "distances=[Borgefors (3,4,5)] output=[16 bits] normalize dynamic=2 connectivity=6");
	run("Set Label Map", "colormap=Spectrum background=Black shuffle");
	rename("temp");
	
	//Size filter and create binary
	run("Label Size Filtering", "operation=Greater_Than size=" + filter_nuclei);
	setThreshold(1.0000, 1000000000000000000000000000000.0000);
	run("Convert to Mask", "background=Dark black");
	
	
	rename("nuclei");
	
	calibrate_image(raw_image, "nuclei");
	
	//Close
	close("nuclei_raw");
	close("temp");
}


//Segment spot marker
function segment_spots_cells(raw_image,spot_channel) { 
	//Isolate channel
	selectWindow(raw_image);
	run("Duplicate...", "duplicate channels=" + spot_channel + "-" + spot_channel);
	rename("spot_raw");
	
	//Normalize the entire stack
	run("Enhance Contrast...", "saturated=0 normalize process_all use"); 
	
	// Init GPU
	run("CLIJ2 Macro Extensions", "cl_device=[NVIDIA GeForce RTX 3090]");
	Ext.CLIJ2_clear();
	
	// difference of gaussian
	image_1 = "spot_raw";
	Ext.CLIJ2_pushCurrentZStack(image_1);
	
	// Difference Of Gaussian2D
	sigma1x = 2;
	sigma1y = 2;
	sigma2x = 10;
	sigma2y = 10;
	Ext.CLIJ2_differenceOfGaussian2D(image_1, image_2, sigma1x, sigma1y, sigma2x, sigma2y);
	
	// pull dif of gaussian image
	Ext.CLIJ2_pull(image_2);
	
	// Cleanup by the end
	Ext.CLIJ2_clear();
	
	//threshold
	resetMinAndMax();
	setAutoThreshold("Moments dark no-reset stack");
	run("Convert to Mask", "method=Moments background=Dark black");
	rename("spots");
	
	calibrate_image(raw_image, "spots");
	
	//Segmenting all cell signal using the spot marker
	selectWindow("spot_raw");
	run("Gaussian Blur...", "sigma=10 stack");
	setAutoThreshold("MinError dark no-reset stack");
	run("Convert to Mask", "method=MinError background=Dark black");
	run("Fill Holes", "stack");
	rename("totalcell_binary");
}

//Segment cluster marker
function segment_cluster(raw_image,cluster_channel,threshold) { 
	//Isolate channel
	selectWindow(raw_image);
	run("Duplicate...", "duplicate channels=" + cluster_channel + "-" + cluster_channel); // In the case of Jacqueline 
	rename("cluster_raw");
	
	//Normalize the entire stack
	run("Enhance Contrast...", "saturated=0 normalize process_all use"); 
	
	//filter and background substraction
	run("Median...", "radius=5 stack");
	run("Subtract Background...", "rolling=50 stack");
	
	//threshold
	setThreshold(threshold_cluster, 65535, "raw");
	run("Convert to Mask", "method=Otsu background=Dark black");
	
	rename("cluster");
}

// Add an additional column with the metadata
//It takes as an input a table, a new column name and a value and adds
// a new column with the string repeated. It is usually useful for metadata.
function add_metadata(table_name,newcolumn_name,value) {
	selectWindow(table_name);
	rows = getValue("results.count");
	newcolumn_array = newArray(rows);
	for (row = 0; row < rows; row++) {
		newcolumn_array[row] = value;
	}
	Table.setColumn(newcolumn_name, newcolumn_array);
}

//Get the histogram of a stack with the counts in fractions of 1 additionally
function getstackhisto(image, binsize, pixelcount,id) {
	setOption("ExpandableArrays", true);
	//Select image of interest
	selectWindow(image);
	
	//Transform the stack to a 2D image
	run("Make Montage...", "columns=1 rows=" + nSlices + " scale=1");
	
	//Calculate the different values to feed into the histogram function
	selectWindow("Montage");
	resetMinAndMax();
	getMinAndMax(min, max);
	print(max);
	//max = (binsize * Math.ceil(max / binsize)) + binsize; // Round up max to the nearest multiple of binsize + binsize to round up the range
	max = binsize * Math.ceil(max / binsize);
	
	bins = max / binsize;
	
	// Initialize the resulting array and the counts_bis array
	result_fraction = newArray(bins);
	result_count = newArray(bins);
	id_array = newArray(bins);
	
	// Process montage
	getHistogram(values, counts, bins, 0, max);
	
	//Get the fraction array
	//Control for the fact that some histogram calculations take the background 0 as signal too.
	sum_fraction = 0;
	for (i = 0; i < values.length; i++) {
		result_fraction[i] = counts[i]/pixelcount;
		id_array[i] = id;
		if (i > 0) {
			sum_fraction = sum_fraction + result_fraction[i]; 
		}
	}
	
	result_fraction[0] = 1-sum_fraction;
	counts[0] = result_fraction[0] * pixelcount;
	
	//Generate a results table containing the stack histogram
	Table.create("Results");
	Table.setColumn("Values", values);
	Table.setColumn("Counts", counts);
	Table.setColumn("Counts_fraction", result_fraction);
	add_metadata("Results","ROI_ID",id);
	
	close("Montage");
}

//Generate distance image
function isolate(image,roi) { 
	selectWindow(image);
	roiManager("Select", roi);
	run("Duplicate...", "duplicate");
	rename(roi);
	run("Clear Outside", "stack");
}

// This function accepts a 32-bitNND image to be scaled by the min max of distance_image
function distance_scaling(image_name,distance_image) { 
	selectWindow(distance_image);
	
	//Transform the stack to a 2D image
	run("Make Montage...", "columns=1 rows=" + nSlices + " scale=1");
	resetMinAndMax();
	rename("test");
	
	getMinAndMax(min, max);
	distance_delta = round(max-min);
	
	selectWindow(image_name);
	run("Subtract...", "value=" + min + " stack");
	run("Divide...", "value=" + distance_delta + " stack");
	close("test");
}

// This function takes a roi list with annotated cells. Isolates the totalcell_binary, spots and cluster per ROI. Calculates the local separation of the spots image and the NND of the spots binary to the cluster and generates a separate table per. 
function process_roi() { 
	setOption("ExpandableArrays", true);
	roi_count = roiManager("count");
	for (roi = 0; roi < roi_count; roi++) {
	//for (roi = 0; roi < test_roi; roi++) {
		//Isolate single cells
		isolate("totalcell_binary",roi);
		rename("cell_" + roi);
		
		//Isolate single nuclei
		isolate("nuclei",roi);
		rename("nucleus_" + roi);
		
		//Generate single cell cytosol
		imageCalculator("Subtract create stack", "cell_" + roi, "nucleus_" + roi);
		rename("cytosol_" + roi);
		
		//Get voxel count for detecting cells with and without clusters
		run("Analyze Regions 3D", "voxel_count surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");		
		cell_voxelcount = Table.get("VoxelCount", 0);
		
		//Isolate single cell spots
		isolate("spots",roi);
		rename("spots_" + roi);
		
		//Get voxel count for fraction normalization in spatial and morphometric measurements
		run("Analyze Regions 3D", "voxel_count surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
		spots_voxelcount = Table.get("VoxelCount", 0);
		
		//Isolate single cell cluster
		isolate("cluster",roi);
		rename("cluster_" + roi);
		
		//Get voxel count for detecting cells with and without clusters
		run("Analyze Regions 3D", "voxel_count surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");		
		Table.rename("cluster_" + roi + "-morpho", "Results");
		
		cluster_voxelcount = 0;
		if (nResults > 0) {
			cluster_voxelcount = Table.get("VoxelCount", 0);
		}

		//Generate the inverted image for the separation calculation
		imageCalculator("Subtract create stack", "cytosol_" + roi,"spots_" + roi);
		rename("spot_inverted_" + roi);
		
		//Get voxel count for fraction normalization in spatial and morphometric measurements
		run("Analyze Regions 3D", "voxel_count surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
		spot_inverted_voxelcount = Table.get("VoxelCount", 0);
		
		//Local thickness of the separation binary
		run("Local Thickness (masked, calibrated, silent)");
		rename("spot_localsep" + roi);
		calibrate_image(raw_image, "spot_localsep" + roi);
		
		//Get distance from cluster
		calibrate_image(raw_image, "cluster_" + roi);
		selectWindow("cluster_" + roi);
		run("Distance Transform 3D");
		rename("cluster_distance" + roi);
		
		//Generate spot distance to cluster image
		selectWindow("spots_" + roi);
		run("Replace/Remove Label(s)", "label(s)=255 final=1");//transform the binary to label 1 for multiplication
		rename("spots_255to1_" + roi);
		imageCalculator("Multiply create 32-bit stack", "spots_255to1_" + roi,"cluster_distance" + roi);//Generate the NND to cluster from spots image by multiplying the spots label 1Xcluster distance image
		rename("spot_nndtocluster" + roi);
		
		//Get the histogram counts and fraction of local separation
		getstackhisto("spot_localsep" + roi, binsize_sep, spot_inverted_voxelcount,roi);
		add_metadata("Results","Aggregate_ratiotoCell",cluster_voxelcount/cell_voxelcount);
		add_metadata("Results","Filename",name);
		
		//Save results
		saveAs("Results", dir + "/Analysis/" + name + "_" + roi + "_localsep.csv");
		
		//Get the histogram counts and fraction of distance to cluster
		getstackhisto("spot_nndtocluster" + roi, binsize_dist, spots_voxelcount,roi);
		add_metadata("Results","Aggregate_ratiotoCell",cluster_voxelcount/cell_voxelcount);
		add_metadata("Results","Filename",name);
		
		//Save results
		saveAs("Results", dir + "/Analysis/" + name + "_" + roi + "_nndtocluster.csv");
		
		//Scale the distance to make cells of different sizes comparable
		selectImage("cytosol_" + roi);
		run("Subtract...", "value=254 stack");
		
		imageCalculator("Multiply create 32-bit stack", "cytosol_" + roi,"cluster_distance" + roi);//Generate the NND to cluster from cell
		rename("cytosol_nndtocluster" + roi);
		
		//
		distance_scaling("spot_nndtocluster" + roi,"cytosol_nndtocluster" + roi);
		
		//Get the histogram counts and fraction of distance to cluster
		getstackhisto("spot_nndtocluster" + roi, binsize_dist_scaled, spots_voxelcount,roi);
		add_metadata("Results","Aggregate_ratiotoCell",cluster_voxelcount/cell_voxelcount);
		add_metadata("Results","Filename",name);
	
		//Save results
		saveAs("Results", dir + "/Analysis/" + name + "_" + roi + "_nndtocluster_scaled.csv");
		
		//Close opened images
		/**/
		close("cell_" + roi);
		close("nucleus_" + roi);
		close("cytosol_" + roi);
		close("spots_" + roi);
		close("spot_inverted_" + roi);
		close("spot_localsep" + roi);
		close("spots_255to1_" + roi);
		close("spot_nndtocluster" + roi);
		close("cytosol_nndtocluster" + roi);
		close("cluster_" + roi);
		close("cluster_distance" + roi);
		
		//close("\\Others");
		
		close("cytosol_" + roi + "-morpho");
		close("cluster_" + roi + "-morpho");
		close("spots_" + roi + "-morpho");
		close("spot_inverted_" + roi + "-morpho");
	}
}
