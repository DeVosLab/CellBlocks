/*	CellBlocks_v22.ijm
	***********************
	Author: 			Winnok H. De Vos
	Modified by: 		Winnok H. De Vos & Marlies Verschuuren
	Date Created: 		November 23rd, 2007
	Date Last Modified:	December 3th, 2024 
 	
 	Description: 
 	Macro Set dedicated to cytometric analysis of n-dimensional images containing at least one channel 
 	of a nuclear dye to be used for segmentation of individual cells.
 	It returns nuclear as well as cellular shape descriptors, intensity and spot parameters of compartments 
 	in selected channels. The cytoplasm is segmented blindly, i.e. from the nuclear rois by region growing or using a cellular counterstain.
	N-dimensional files should be saved as separate regions from a mosaic acquisition. If they are not separated, first use the separate[r] command.	
	Requires FeatureJ and imagescience plugins (E. Meijering): download at http://www.imagescience.org/meijering/software/featurej/
	From v6 onwards: also requires Stardist and CSBDeep plugins (overlapping nuclei segmentation): https://imagej.net/StarDist
	From v11 onwards: also requires GLCM2 plugin for texture analysis: download at https://github.com/miura/GLCM2/blob/master/target/GLCM2_-1.0.0.jar
	From v15 onwards: also requires Trained_Cell_Segmentation.bsh plugin and a trained model (*.model)
	From v16 onwards: also requires PTBIOP plugin for cellpose and a working cellpose environment: https://github.com/BIOP/ijl-utilities-wrappers and https://github.com/MouseLand/cellpose#local-installation
	From v16 onwards: also requires Cellpose_Script_Wrapper.groovy plugin for cellpose
	
	Change Log
	+ Included a simple method to determine the distance from the nucleus
	+ Included Fixed threshold settings for nuclei and spots
	+ Added dilation for cytoplasm measurements as alternative regions of influence
	+ also incldued an option to just measure cellular regions without nuclear area
	+ introduced laplacian detector to detect nuclei in very noisy images
	+ included cell segmentation for a cellular counterstain
	+ resolved discrepancy between nuc and cell segmentation
	+ included fold analysis for lamin channel
	+ added a border removal (rim) and flatfield correction
	+ added a bidimensional plot (two params) with gating function to select cells based on cutoff values (only for single images) (150918)
	+ added interactive lookup/selection tool to check ROI with plot point and vice versa (170918)
	+ corrected an error that caused cell roi detection to be incorrect when excedding 255 objects (changed the pid to be 16-bit instead of 8-bit) (170918)
	+ corrected an error in single image and batch mode that caused bug upon lack of spots (041119)
	+ added a nuclear shrink option (to allow internal rim)(041119)
	+ included Stardist for improved nuclei detection (trained network) and removed nuclear stringency check (040620)
	+ added simple spot detection without laplace for bulky objects (191020)
	+ added distance calculator between objects(191020)
	+ renamed all variables to facilitate interpretation (191020)
	+ Added a new action tool to separate very large images into tiles of defined size for more convenient analyses (171120)
	+ corrected distance calculation integration and threshold setting for spots with Gauss detection (170321)
	+ restructured the setup and settings display (130421)
	+ added a plaque segmentation module (130402)
	+ corrected fixed threshold setting for plaque segmentation (280621)
	+ corrected retention of plaque results when plaques are absent (250821)
	+ added an additional loop to separate and measure nuclear plaque constituents (260821)
	+ texture analysis (GLCM - haralick) (121021 - Marlies Verschuuren)
	+ add calibration before analysisRegions (121021 - Marlies Verschuuren)
	+ add calibration in mask Stardist analysis (121021 - Marlies Verschuuren)
	+ added edge detection info to nuclear measurement output (261021)
	+ corrected a bug in Gauss spot segmentation method (wrong threshold approach) (261021)
	+ added a fold skeletonization parameter (191121)
	+ added watershed for plaque separation (291121)
	+ added a routine for fusing two channels (in casu DAPI an Calcein) to obtain a better cell channel for segmentation purposes (301121)
	+ extended cell and nucleus measurements to include kurtosis, skewness, shape (cell) (031221)
	+ extended distance measurement between objects from only 'to' to also allow 'within' (e.g. spots inside nucleus) (180222)
 	+ fixed an error that would cause multi-scale segmentation to always use scale 3 (150422)
 	+ corrected stardist normalisation range from 1-100% to 1-99% (260422)
 	+ corrected cell segmentation to correctly assign cytoplasm with nuclei in all cases (else multiple pieces of cytoplasm could be assigned to same nucleus leading to erronous cell-nucleus matching) (131022)
 	+ added applying a model obtained with Trainable Segmentation for cell detection (141222)
	+ changed back the cell separation to initiate from a filtered nuclei mask (instead of full nuclei mask) to avoid oversegmentation (141222)	
	+ cellpose cytoplasm detection (190123 - Marlies Verschuuren)
	+ corrected misassignment (introduced ROI matching) of cells and nuclei when lack of overlap and corrected ROI count issue (300123)
	+ included nuclear channel in cellpose segmentation (080223)
	+ corrected roi-names after ROI matching (080223)
	+ converted mask of positive nuclei to a 16-bit image to allow labelmap of > 255 nuclei (080223)
	+ select none before texture analysis (090323 - Marlies Verschuuren)
	+ fixed an error in split tiles: image size is now independent of portrait/landscape mode (090323 - Marlies Verschuuren)
	+ fixed an error that occurred when using dilation for cell detection: rename Cell_Roi (090323 - Marlies Verschuuren)
	+ fixed an error that occurred using cellpose due to very small pixel regions in the mask (210423 - Marlies Verschuuren)
	+ removed ROI matching and removed check for cell duplicates (v19 - 140623 - Marlies Verschuuren)
	+ add nuclei roi index in cell results and link in summary results (v19 - 140623 - Marlies Verschuuren)
	+ adapt spot index: "exclude_nuclei" + remove "exclude_nuclei" option in summary (v19 - 190623 - Marlies Verschuuren)
	+ work with labeled centroid map of nuclei to link cells (needed for cellpose, since cell borders do not match with voronoi borders) (v19 - 210623 - Marlies Verschuuren)
 	+ fixed error dilation cell (v19b - 240723 - Marlies Verschuuren)
 	+ harmonzied definition of fixed threshold value for nuclei and spot detection (fill in negative max value when using laplace, positive min value when using gauss,  just as displayed for the autothreshold) (v20 - 090724)
 	+ rewrite exclude nucleus section: merge nuc+cytoplasm=cell -> sort roiManager -> XOR if cell>nuc  (v21 - 220724 - Marlies Verschuuren)
 	+ harmonized definition of fixed threshold applied to multi-scale spot detection (v21b - 031024 - Marlies Verschuuren)
 	+ add "microns" identifier to check for calibration (v22 - 031224 - Marlies Verschuuren)
 	+ fix bug spotdetection "indices" not found (v22 - 111224 - Marlies Verschuuren)
 	_________________________________________________________________
*/

/*
 	***********************

	Variable initiation

	***********************
*/


//	String variables
var cells_results					= "";										//	cell region results
var cells_roi_set 					= "";										//	cell region ROIsets
var cells_segmentation_method		= "Threshold";								//  Methods used for cell detection
var cells_threshold					= "Fixed";									//	threshold method for segmentation of cells 
var dir								= "";										//	directory
var distance_results				= "";										//	distance between ROI objects results	
var folds_results					= "";										//	fold results
var folds_roi_set					= "";										//	folds ROI sets 	
var folds_segmentation_method		= "Sobel";									// 	method for segmenting nuclear folds
var folds_threshold					= "Default";								//	threshold method for fold segmentation
var ind_plaques_results				= "";										//	individual plaque results
var log_path						= "";										//	path for the log file
var micron							= getInfo("micrometer.abbreviation");		// 	micro symbol
var model 							= "/Applications/Fiji.app/models/classifier.model"; // location of a trained model for cell segmentation
var nuclei_roi_set 					= "";										//	nuclear ROIsets
var nuclei_results 					= "";										//	nuclear results
var nuclei_segmentation_method		= "Stardist";								// 	Method used for nuclei detection
var nuclei_threshold				= "Default";								//	threshold method for segmentation of nuclei 
var plaques_roi_set 				= "";										//	plaque ROIsets name
var plaques_threshold 				= "Fixed";									//	threshold method for segmentation of plaques 
var plaques_results					= "";										//	plaque results
var order							= "xyczt(default)";							//	hyperstack dimension order
var output_dir						= "";										//	dir for analysis output
var primary_objects					= "Plaques";								//	primary objects to which secondary objects are located (distance calculation)
var results							= "";										//	summarized results	
var secondary_objects				= "Nuclei";									//	Secondary objects for distance calculation
var spots_a_results					= "";										//	spot results ch A
var spots_a_roi_set					= "";										//	spot ROI sets ch A			
var spots_a_segmentation_method		= "Laplace";								//	spots segmentation method ch A
var spots_a_threshold				= "Fixed";									//	threshold method for spot segmentation ch A
var spots_b_results					= "";										//	spot results ch B
var spots_b_roi_set					= "";										//	spot ROI sets ch B	
var spots_b_segmentation_method		= "Laplace";								//	spots segmentation method ch B
var spots_b_threshold				= "Fixed";									//	threshold method for spot segmentation ch B
var suffix							= ".nd2";									//	suffix for specifying the file type

//	Number variables
var channels						= 3;										//	number of channels
var cells_channel					= 3;										//	channel used for segmentation of cells 
var cells_diameter					= 100;										//  approximate diameter of cells 
var cells_filter_scale				= 1;											//	radius for cell smoothing
var cells_fixed_threshold_value		= 1000;										//	fixed maximum threshold for cell segmentation (if auto doesn't work well);
var fields							= 4;										//	number of xy positions
var folds_channel					= 1;										//	channel used for segmentation of folds
var folds_filter_scale				= 1;										//	scale for Laplacian folds
var folds_fixed_threshold_value		= 50;										//	fixed maximum threshold for fold segmentation (if auto doesn't work well);
var folds_min_area					= 2;										//  min fold size
var image_height					= 1000;										//	image height
var image_width						= 1000;										//	image width
var iterations						= 25;										// 	iterations of dilations in region growing (for determining cell boundaries)
var slices							= 7;										//	number of z-slices
var nuclei_channel					= 1;										//	channel used for segmentation of nuclei 
var nuclei_filter_scale				= 0;										// 	radius for nuclei smoothing/laplace
var nuclei_fixed_threshold_value	= 100;										//	fixed maximum threshold for nuclei segmentation (if auto doesn't work well);
var nuclei_min_circularity			= 0.0;										//	min circularity
var nuclei_min_area					= 50;										//	calibrated min nuclear size (in µm2)
var nuclei_max_area					= 5000;										//	calibrated max nuclear size (in µm2)
var nuclei_overlap 					= 0.3;										//	nuclei_overlap amount tolerated for Stardist nuclei detection
var nuclei_probability				= 0.15;										//	minimal nuclei_probability for Stardist nuclei detection
var pixel_size						= 0.3641552;									//	pixel size (µm)
var plaques_channel					= 2;										//	channel used for segmentation of plaques 
var plaques_filter_scale 			= 3;										//	scale for Gauss filtering of plaque channel
var plaques_fixed_threshold_value 	= 200;										//	fixed maximum threshold for plaque segmentation 
var plaques_max_area				= 999999;									//	max plaque size in pixels
var plaques_min_area				= 10;										//	min plaque size in pixels 
var rim								= 0;										//	percentage of the field to analyze
var shrink							= 2;										// 	shrink factor for nuclear contour in fold segmentation
var spots_a_channel					= 1;										//	channel A used for segmentation of spots
var spots_a_filter_scale			= 1;										//	scale for Laplacian spots ch A
var spots_a_fixed_threshold_value	= 180;										//	fixed maximum threshold for spot segmentation 
var spots_a_max_area				= 20;										//	max spot size in pixels ch A
var spots_a_min_area				= 2;										//	min spot size in pixels ch A
var spots_b_channel					= 2;										//	channel B used for segmentation of spots
var spots_b_filter_scale			= 1;										//	scale for Laplacian spots ch B
var spots_b_fixed_threshold_value	= 180;										//	fixed maximum threshold for spot segmentation
var spots_b_max_area				= 20;										//	max spot size in pixels ch B
var spots_b_min_area				= 2;										//	min spot size in pixels ch B
var tile_size						= 2000;										//	size of tiles to separate large mosaic files in

//	Boolean Parameters
var colocalize_spots				= false;									//	colocalize spot ROIs
var distance_within					= true;										//	measure distance of secondary objects within the primary objects (otherwise to)
var exclude_nuclei					= true;										//	analyze cellular regions without nuclear area
var flatfield						= false;									//	perform flatfield correction	
var measure_distances				= false;									//	calculate distances between objects
var nuclei_background				= true;										//	subtract nuclei_background for nuclei segmentation
var nuclei_clahe					= false;									// 	local contrast enhancement for nuclei segmentation
var nuclei_watershed				= false;									//	use nuclei_watershed for nuclei segmentation
var plaques_background				= true;										//	subtract nuclei_background for plaque segmentation
var plaques_clahe					= false;									// 	local contrast enhancement for plaque segmentation
var plaques_watershed				= true;										// 	watershed separation for plaque segmentation
var segment_cells					= false;									//	analyze cell and cytoplasmic ROIs
var segment_folds					= false;									//	analyze fold ROIs
var segment_nuclei					= true;										//	analyze nuclear ROIs (currently redundant)
var segment_plaques					= false;									//	analyze plaques 
var segment_spots					= false;									//	analyze spot ROIs
var texture_analysis 				= false;									//  implement texture analysis
var z_project						= true;										//	generation of max projection images
var skeletonize_folds				= true:										//	skeletonize detected folds to obtain only the thinned outlines

//	Arrays
var cells_segmentation_methods		= newArray("Threshold","Dilation","Trained Model","Cellpose");
var cols 							= newArray("01","02","03","04","05","06","07","08","09","10","11","12");
var dimensions						= newArray("xyczt(default)","xyctz","xytcz","xytzc","xyztc","xyzct");		
var file_list						= newArray(0);								
var file_types 						= newArray(".tif",".tiff",".nd2",".ids",".jpg",".mvd2",".czi");		
var folds_segmentation_methods		= newArray("Gauss","Laplace","Sobel");
var nuclei_segmentation_methods		= newArray("Gauss","Laplace","Stardist");
var objects 						= newArray("Nuclei","Cells","Spots_a","Spots_b","Folds","Plaques");
var prefixes 						= newArray(0);
var rows 							= newArray("A","B","C","D","E","F","G","H");
var spot_segmentation_methods		= newArray("Gauss","Laplace","Multi-Scale");
var threshold_methods				= getList("threshold.methods");	
var threshold_methods				= Array.concat(threshold_methods,"Fixed");	

/*
 	***********************

		Macros

	***********************
*/

macro "Autorun"
{
	erase(1);
}
/*
macro "Split Stack Into Fields Action Tool - C888 R0077 R9077 R9977 R0977"
{
	setBatchMode(true);
	splitRegions();
	setBatchMode("exit and display");
}

macro "Split Large File Into Tiles Action Tool - Cf88 R0077 C888 R9077 R9977 R0977"
{
	setBatchMode(true);
	splitTiles();
	setBatchMode("exit and display");
}
*/

macro "Z Project Action Tool - C888 R0099 R3399 R6699 R9999"
{
	setBatchMode(true);
	projectImages();
	setBatchMode("exit and display");
}
/*
macro "Sharpest Z Selection Action Tool - C888 R0099 R3399 Cf88 R6699 C888 R9999"
{
	setBatchMode(true);
	selectImages();
	setBatchMode("exit and display");
}

macro "Add Fused Cell Channel Action Tool - C888 R0099 R3399 Cf88 R6699 C888 R9999"
{
	setBatchMode(true);
	fuseImages();
	setBatchMode("exit and display");
}
*/

macro "Setup Action Tool - C888 T5f16S"
{
	setup();
}

macro "Analyse Single Image Action Tool - C888 T5f161"
{
	erase(0);
	setBatchMode(true);
	dir = getInfo("image.directory");
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir))File.makeDirectory(output_dir);
	start = getTime();
	title = getTitle; 
	prefix = substring(title,0,lastIndexOf(title,suffix));
	setFileNames(prefix);
	id = getImageID;
	if(flatfield)id = flatfield_correct(id);
	mid = segmentNuclei(id,nuclei_channel,1); 
	nuclei_nr = roiManager("count");
	if(nuclei_nr>0)roiManager("Save",nuclei_roi_set);
	if(nuclei_nr>0 && segment_cells)
	{
		cell_nr = segmentRegions(id, mid, cells_channel, iterations);
		if(cell_nr>0)roiManager("Save",cells_roi_set);
		else {File.delete(nuclei_roi_set); nuclei_nr=0;} 
		roiManager("reset");
	}
	if(isOpen(mid)){selectImage(mid); close;}
	if(nuclei_nr>0 && segment_folds)
	{
		if(roiManager("count")==0)roiManager("Open",nuclei_roi_set);
		foldNr = segmentFolds(id, folds_channel);
		if(foldNr>0)
		{
			roiManager("Save",folds_roi_set);
			roiManager("reset");
		}
	}
	if(nuclei_nr>0 && segment_spots)
	{
		roiManager("reset");
		if(spots_a_channel>0)
		{
			args	= newArray(spots_a_channel,spots_a_segmentation_method,spots_a_threshold,spots_a_fixed_threshold_value,spots_a_filter_scale,spots_a_min_area,spots_a_max_area);
			snr 	= segmentSpots(id,spots_a_channel,args);
			if(snr>0)
			{
				roiManager("Save",spots_a_roi_set);
				roiManager("reset");
			}
		}
		if(spots_b_channel>0)
		{
			args	= newArray(spots_b_channel,spots_b_segmentation_method,spots_b_threshold,spots_b_fixed_threshold_value,spots_b_filter_scale,spots_b_min_area,spots_b_max_area);
			snr 	= segmentSpots(id,spots_b_channel,args);
			if(snr>0)
			{
				roiManager("Save",spots_b_roi_set);
				roiManager("reset");
			}
		}
	}
	if(nuclei_nr>0 && segment_plaques)
	{
		roiManager("reset");
		plaque_density = 0;
		plaque_density = segmentPlaques(id, plaques_channel);
		if(plaque_density>0)
		{
			roiManager("Save",plaques_roi_set);
			roiManager("reset");
		}
	}
	readout = analyzeRegions(id);
	if(readout)summarizeResults();
	print((getTime()-start)/1000,"sec");
	print("Analysis Done");
	setBatchMode("exit and display");
}

macro "Batch Analysis Action Tool - C888 T5f16#"
{
	erase(1);
	setBatchMode(true);
	setDirectory();
	prefixes = scanFiles();
	fields = prefixes.length;
	setup();
	start = getTime();
	for(i=0;i<fields;i++)
	{
		prefix = prefixes[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id 		= getImageID;
		if(flatfield)id = flatfield_correct(id);
		mid 	= segmentNuclei(id,nuclei_channel,1); 
		nuclei_nr 	= roiManager("count");
		if(nuclei_nr>0)roiManager("Save",nuclei_roi_set);
		if(nuclei_nr>0 && segment_cells)
		{
			cell_nr = segmentRegions(id, mid, cells_channel, iterations);
			if(cell_nr>0)roiManager("Save",cells_roi_set);
			else {File.delete(nuclei_roi_set); nuclei_nr=0;} 
			roiManager("reset");
		}
		if(isOpen(mid)){selectImage(mid); close;}
		if(nuclei_nr>0 && segment_folds)
		{
			if(roiManager("count")==0)roiManager("Open",nuclei_roi_set);
			foldNr = segmentFolds(id, folds_channel);
			if(foldNr>0)
			{
				roiManager("Save",folds_roi_set);
				roiManager("reset");
			}
		}
		if(nuclei_nr>0 && segment_spots)
		{
			roiManager("reset");
			if(spots_a_channel>0)
			{
				args	= newArray(spots_a_channel,spots_a_segmentation_method,spots_a_threshold,spots_a_fixed_threshold_value,spots_a_filter_scale,spots_a_min_area,spots_a_max_area);
				snr 	= segmentSpots(id,spots_a_channel,args);
				if(snr>0)
				{
					roiManager("Save",spots_a_roi_set);
					roiManager("reset");
				}
			}
			if(spots_b_channel>0)
			{
				args	= newArray(spots_b_channel,spots_b_segmentation_method,spots_b_threshold,spots_b_fixed_threshold_value,spots_b_filter_scale,spots_b_min_area,spots_b_max_area);
				snr 	= segmentSpots(id,spots_b_channel,args);
				if(snr>0)
				{
					roiManager("Save",spots_b_roi_set);
					roiManager("reset");
				}
			}
		}
		if(nuclei_nr>0 && segment_plaques)
		{
			roiManager("reset");
			plaque_density = 0;
			plaque_density = segmentPlaques(id, plaques_channel);
			if(plaque_density>0)
			{
				roiManager("Save",plaques_roi_set);
				roiManager("reset");
			}
		}
		readout = analyzeRegions(id);
		if(readout)summarizeResults();
		selectImage(id); close;
		erase(0);
	}
	print((getTime()-start)/1000,"sec");
	if(isOpen("Log")){selectWindow("Log");saveAs("txt",log_path);}
	print("Complete Analysis Done");
	setBatchMode("exit and display");
}


macro "Segment Nuclei Action Tool - C999 H11f5f8cf3f181100 C999 P11f5f8cf3f18110 Ceee V5558"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	c 		= getNumber("Nuclear Channel",nuclei_channel);
	mid 	= segmentNuclei(id,c,0);
	setBatchMode("exit and display");
}

macro "Segment Cells Action Tool - C999 H11f5f8cf3f181100 Ceee P11f5f8cf3f18110"
{
	erase(0);
	setBatchMode(true);
	id = getImageID;
	c1 = getNumber("Nuclear Channel",nuclei_channel);
	c2 = getNumber("Cellular Mask Channel",cells_channel); // 0 if no additional mask
	mid = segmentNuclei(id,c1,1); 
	nuclei_nr = roiManager("count");
	if(nuclei_nr>0)cell_nr = segmentRegions(id, mid, c2, iterations);
	setBatchMode("exit and display");
}

macro "Segment Folds Action Tool - C999 H11f5f8cf3f181100 C999 P11f5f8cf3f18110 Ceee V3329 V6519 V8519 Va827"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	c1 		= getNumber("Nuclear Channel",nuclei_channel);
	c2 		= getNumber("Fold Channel",folds_channel); 
	mid 	= segmentNuclei(id,c1,0);
	fnr 	= segmentFolds(id,c2);
	setBatchMode("exit and display");
}

macro "Segment Spots Action Tool - C999 H11f5f8cf3f181100 C999 P11f5f8cf3f18110 Ceee V3633 V4b33 V7633 Va933"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	c 		= getNumber("Spot Channel",spots_a_channel);
	if(c == spots_a_channel){args = newArray(spots_a_channel,spots_a_segmentation_method,spots_a_threshold,spots_a_fixed_threshold_value,spots_a_filter_scale,spots_a_min_area,spots_a_max_area);}
	if(c == spots_b_channel){args = newArray(spots_b_channel,spots_b_segmentation_method,spots_b_threshold,spots_b_fixed_threshold_value,spots_b_filter_scale,spots_b_min_area,spots_b_max_area);}
	snr 	= segmentSpots(id,c,args);
	setBatchMode("exit and display");
}

macro "Segment Plaques Action Tool - C999 R00FF F00ff Ceee V1158 V9758"
{
	erase(0);
	setBatchMode(true);
	id = getImageID;
	c = getNumber("Plaque Channel",plaques_channel); 
	plaque_density = segmentPlaques(id, c);
	setBatchMode("exit and display");
}

macro "Toggle Overlay Action Tool - Caaa O11ee"
{
	toggleOverlay();
}

macro "[t] Toggle Overlay"
{
	toggleOverlay();
}

macro "Verification Stack Action Tool - C888 T5f16V"
{
	erase(1);
	setBatchMode(true);
	setDirectory();
	prefixes = scanFiles();
	subset = getBoolean(""+prefixes.length+" images, select subset?");
	if(subset)
	{
		//	Compose an array containing all 96 well labels
		allWells = newArray(0);	
		for(i = 0; i < rows.length; i++)
  		{
  			for(j = 0; j < cols.length; j++)
  			{
  				allWells = Array.concat(allWells, rows[i] + cols[j]);
  			}
  		}
  		selection = "";
  		defaults = newArray(allWells.length);
  		Array.fill(defaults,0);		
  		selWells = newArray(allWells.length);		
		Dialog.create("Verification Stack");
		Dialog.addMessage("Select wells to visualize...")
		Dialog.addCheckboxGroup(rows.length,cols.length,allWells,defaults);
  		Dialog.show;
  		for(i = 0; i < rows.length; i++)
  		{
  			for(j = 0; j < cols.length; j++)
  			{
  				v = Dialog.getCheckbox();
  				if(v>0)selection = selection + allWells[i*cols.length+j];
  			}
  		}
		names = newArray(0);
  		for(i = 0;i<prefixes.length;i++)
		{
			name = prefixes[i];
			well = substring(name,indexOf(name,"Well")+4,indexOf(name,"Well")+7);
			if(indexOf(selection,well)>-1)names = Array.concat(names,name);
		}
	}
	else names = prefixes;
	createOverlay(names);
	setBatchMode("exit and display");
	run("Channels Tool... ");
}

macro "Gate Results Action Tool - C888 T5f16G"
{
	if(nImages==0)open();
	id = getImageID; 
	selectImage(id);
	title = getTitle;
	if(roiManager("count")==0)
	{
		dir = getDirectory("image");
		output_dir = dir+"Output"+File.separator;
		prefix = substring(title,0,lastIndexOf(title,suffix));
		nuclei_roi_set = output_dir+prefix+"_nuclei_roi_set.zip";
		results	= output_dir+prefix+"_summary.txt";
		roiManager("Open",nuclei_roi_set);
		run("Results... ", "open=["+results+"]");		
	}
	run("Select None");
	roiManager("deselect");
	// extract result labels
	nr = nResults;
	selectWindow("Results");
	lines = split(getInfo(),'\n');
	labels = split(lines[0],'\t');
	labels = Array.concat(labels,"index");
	// select metrics to gate (index is no gate)
	Dialog.create("Gate by Metric (highlight within selected range)");
	Dialog.addChoice("Metric_X",labels);
	Dialog.addNumber("Lower Value",0);
	Dialog.addNumber("Upper Value","Infinity");
	Dialog.addChoice("Metric_Y",labels);
	Dialog.addNumber("Lower Value",0);
	Dialog.addNumber("Upper Value","Infinity");
	Dialog.show
	metric_X		= Dialog.getChoice();
	lo_level_X		= Dialog.getNumber;
	up_level_X		= Dialog.getNumber;
	metric_Y		= Dialog.getChoice();
	lo_level_Y		= Dialog.getNumber;
	up_level_Y		= Dialog.getNumber;
	// apply gating
	index = 0;
	pos_values_X = newArray(0);
	neg_values_X = newArray(0);
	pos_values_Y = newArray(0);
	neg_values_Y = newArray(0);	
	selectImage(id);
	run("From ROI Manager");
	run("Overlay Options...", "stroke=gray width=1 fill=none apply");
	run("Colors...", "foreground=white nuclei_background=black selection=red");
	if(lo_level_X+lo_level_Y==NaN && up_level_X+up_level_Y==NaN)
	{
		gate = 0; 
		print("no gating");
	}
	else 
	{
		gate = 1;
	}
	min_X = 0; max_X = 0; min_Y = 0; max_Y = 0; //plot limits
	for(i=0;i<nr;i++)
	{
		if(metric_X == "index")value_X = i+1;
		else value_X = getResult(""+metric_X,i);
		min_X = minOf(min_X, value_X); max_X = maxOf(max_X, value_X); 
 		if(metric_Y == "index")value_Y = i+1;
		else value_Y = getResult(""+metric_Y,i);
		min_Y = minOf(min_Y, value_Y); max_Y = maxOf(max_Y, value_Y); 
		if(value_X>=lo_level_X && value_X<=up_level_X && value_Y>=lo_level_Y && value_Y<=up_level_Y)
		{
			pos_values_X = Array.concat(pos_values_X,value_X);
			neg_values_X = Array.concat(neg_values_X,value_X);
			pos_values_Y = Array.concat(pos_values_Y,value_Y);
			neg_values_Y = Array.concat(neg_values_Y,value_Y);
			if(gate)
			{
				setKeyDown("shift");
				roiManager("select",i);
			}
			index++;
		}
		else 
		{
			pos_values_X = Array.concat(pos_values_X,0/0);
			neg_values_X = Array.concat(neg_values_X,value_X);
			pos_values_Y = Array.concat(pos_values_Y,0/0);
			neg_values_Y = Array.concat(neg_values_Y,value_Y);
		}
	}
	// Create scatterplot
	Plot.create("Plot", ""+metric_X,""+metric_Y);
	Plot.setColor("black","lightGray");
	Plot.add("circles", neg_values_X,neg_values_Y);
	Plot.setColor("black","red");
	Plot.add("circles", pos_values_X,pos_values_Y);
	Plot.setColor("red");
	Plot.drawLine(min_X, lo_level_Y, max_X, lo_level_Y);
	Plot.drawLine(min_X, up_level_Y, max_X, up_level_Y);
	Plot.drawLine(lo_level_X, min_Y, lo_level_X, max_Y);
	Plot.drawLine(up_level_X, min_Y, up_level_X, max_Y);
	Plot.setLimitsToFit();
	Plot.show();
	// return numbers
	print(index,"/",nResults,"positive objects");
}

macro "Selection Tool - C888 L8085 L0858 L8b8f Lb8f8 Cf88 O7722"
{
	// interactive tool
  	d = 100000;  	// distance metric
	sel = -1;		//	index of selection roi
	x_values = newArray(0);
	y_values = newArray(0);
	plot2roi = 0;
	getCursorLoc(x, y, z, flags);
    toScaled(x, y);
    if(flags==16)
    {
   		if(getInfo("window.type")=="Plot")plot2roi = 1;
   		if(plot2roi) 
   		// find corresponding ROI in image
   		{
   			Plot.getValues(x_values, y_values);
   		}
  		else 
  		// find corresponding point in scatterplot
  		{
  			for(r = 0;r < nResults; r++)
  			{
  				x_values = Array.concat(x_values,getResult("X",r));
  				y_values = Array.concat(y_values,getResult("Y",r));
  			}
  		}
  		for (i=0; i<x_values.length; i++)
    	{
    		x2 = x_values[i];
			y2 = y_values[i];
			di = sqrt(pow(x-x2,2)+pow(y-y2,2));
			if(di<d)
			{
				d = di;
				sel = i; 
			}
		}
		print("ROI index:",sel); 
		i = 1;
		if(plot2roi)
		{
			while(i<=nImages && getInfo("window.type")=="Plot")
			{
				selectImage(i); 
				i++;
			}
			Overlay.hide;
			roiManager("select",sel);
		}
		else 
		{
			roiManager("select",sel);
			while(i<=nImages && getInfo("window.type")!="Plot")
			{
				selectImage(i); 
				i++;
			}
			Plot.getValues(x_param, y_param);
			x_coord = x_param[sel];  // get the correct coord in the scatterplot
			y_coord = y_param[sel];
			toUnscaled(x_coord,y_coord);
			makePoint(x_coord,y_coord);
			run("Properties... ", "  stroke=cyan point=Cross size=[Extra Large]");
		}		
	}
}

/*
 	***********************

		Functions

	***********************
*/

function splitRegions()
{
	erase(1);
	Dialog.create("Split Fields...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Image",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.addNumber("Channels",3);
	Dialog.addNumber("Frames",5);
	Dialog.addNumber("Slices",1);
	Dialog.addChoice("Dimension order",dimensions, order);
	Dialog.addCheckbox("Z-Project",z_project);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	channels	= Dialog.getNumber;
	fields 		= Dialog.getNumber;
	slices 		= Dialog.getNumber;
	order 		= Dialog.getChoice;
	z_project 	= Dialog.getCheckbox;

	dir = getDirectory("");
	file_list = getFileList(dir);
	destination_dir = dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			run("Stack to Hyperstack...", "order="+order+" channels="+channels+" slices="+slices+" frames="+fields+" display=Color");
			id = getImageID;
			title = getTitle; 
			for(f=1;f<=fields;f++)
			{
				selectImage(id);
				Stack.setFrame(f);
				//sub = getInfo("image.subtitle"); 
				//if(sub==""){print("no subtitle, using image title"); name = title;}
				//else name = substring(sub,indexOf(sub,"- ")+1,lastIndexOf(sub,ext));
				//print(name); 
				name = getTitle;
				run("Duplicate...", "title=copy duplicate frames="+f);
				cid = getImageID;
				selectImage(cid);	
				if(z_project && slices>1)
				{
					selectImage(cid);
					run("Z Project...", " projection=[Max Intensity]");
					//run("Z Project...", " projection=[Standard Deviation]");
					selectImage(cid); close;
					selectWindow("MAX_copy");
					cid = getImageID;
				}
				selectImage(cid);
				if(f<10)nr="0"+f; 
				else nr=f;
				saveAs(suffix,destination_dir+pre+name+"_xy"+nr+suffix);
				close;
			}
			selectImage(id); close;
		}
	}
	print("Done");
}

function splitTiles()
{
	erase(1);
	Dialog.create("Split Tiles...");
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.addNumber("Tile size (pixels)", tile_size);
	Dialog.show;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	tile_size 	= Dialog.getNumber;

	dir = getDirectory("");
	file_list = getFileList(dir);
	destination_dir = dir+"Tiles"+File.separator;
	File.makeDirectory(destination_dir);
	
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			id = 0;
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default view=Hyperstack ");
			id 		= getImageID;
			title 	= getTitle;
			name = substring(title,0,indexOf(title,ext));
			print(i, name);
			im_width = getWidth;
			im_height = getHeight();
			im_size = minOf(im_width, im_height); // Take smallest size
			tile_nr = Math.ceil(im_size/tile_size);
			for(x=0;x<tile_nr;x++)
			{
				for(y=0;y<tile_nr;y++)
				{
					x_start = x*tile_size;
					y_start = y*tile_size;
					selectImage(id);
					run("Specify...", "width="+tile_size+" height="+tile_size+" x="+x_start+" y="+y_start+" slice=1");
					run("Duplicate...","title=copy duplicate");
					saveAs(suffix,destination_dir+name+"_x"+x+"y"+y+suffix);
					close;
				}
			}
			selectImage(id); close;
		}
	}
	print("Done");
}

function projectImages()
{
	erase(1);
	Dialog.create("Project Images...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Image",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	dir 		= getDirectory("");
	file_list 	= getFileList(dir);
	destination_dir 	= dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			print(i+1);
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			ns = nSlices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+ns/channels+" frames=1 display=Color");
			id = getImageID;
			title = getTitle;
			run("Z Project...", "projection=[Max Intensity]");
			zid = getImageID;		
			selectImage(zid); saveAs(suffix,destination_dir+pre+title+suffix);
			selectImage(id); close;
			selectImage(zid); close;
		}
	}
	print("Done");
}

function selectImages()
{
	erase(1);
	Dialog.create("Z-Select Images...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Image",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	dir 		= getDirectory("");
	file_list 		= getFileList(dir);
	destination_dir 	= dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	for(i = 0;i < file_list.length; i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			print(i+1);
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			id 		= getImageID;
			title 	= getTitle;
			getDimensions(width, height, channels, slices, frames);
			run("Duplicate...","title=Div duplicate");
			did 	= getImageID;
			run("Find Edges", "stack");
			ss = newArray(channels);	//	sharpest slices
			for(c = 1; c <= channels; c++)
			{
				stdmax = 0;
				selectImage(did);
				Stack.setChannel(c);
				for(z = 1; z <= slices; z++)
				{
					Stack.setSlice(z);
					getRawStatistics(nPixels, mean, min, max, std);
					if(std>stdmax)
					{
						stdmax = std;
						ss[c-1] = z;
					}
				}
			}
			selectImage(did); close;		
			for(c = 1; c <= channels; c++)
			{
				selectImage(id);
				Stack.setChannel(c);
				Stack.setSlice(ss[c-1]);
				if(c==1)
				{
					run("Duplicate...","title=Sel duplicate channels="+c+" slices="+ss[c-1]);
					zid = getImageID;
				}
				else
				{
					run("Select All");
					run("Copy");
					selectImage(zid);
					run("Add Slice");
					run("Paste");
					run("Select None");
				}
			}					
			selectImage(zid); saveAs(suffix,destination_dir+pre+title+suffix);
			selectImage(id); close;
			selectImage(zid); close;
		}
	}
	print("Done");
}

function fuseImages()
{
	erase(1);
	Dialog.create("Add Fake Cell Channel To Images by Fusing 2 Channels...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Extended",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.addNumber("Channel 1",1,1,0,"");
	Dialog.addNumber("Channel 2",2,1,0,"");
	Dialog.show;
	dest 				= Dialog.getString;
	pre					= Dialog.getString;
	ext					= Dialog.getChoice;
	suffix 				= Dialog.getChoice;
	channel_1			= Dialog.getNumber;
	channel_2			= Dialog.getNumber;	
	dir 				= getDirectory("");
	file_list 			= getFileList(dir);
	destination_dir 	= dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			print(i+1);
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			id = getImageID;
			selectImage(id);
			title = getTitle; 
			run("Duplicate...", "title=copy duplicate channels="+channel_1+"-"+channel_2);
			copy_id = getImageID;
			selectImage(copy_id);
			run("Properties...", "channels=1 slices=2 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
			run("Z Project...", "projection=[Sum Slices]"); 
			selectImage(copy_id); close;
			selectWindow("SUM_copy");
			run("Select All");
			run("Copy");
			close; 
			selectImage(id);
			getDimensions(width, height, channels, slices, frames);
			Stack.setChannel(channels);
			run("Add Slice", "add=channel");
			Stack.setChannel(channels+1);
			run("Paste");
			selectImage(id); saveAs(suffix,destination_dir+pre+title+suffix);
			close;
		}
	}
	print("Done");
}

function setOptions()
{
	run("Options...", "iterations=1 count=1");
	run("Colors...", "foreground=white nuclei_background=black selection=yellow");
	run("Overlay Options...", "stroke=red width=1 fill=none");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
}

function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}

function erase(all)
{
	if(all){
		print("\\Clear");
		run("Close All");
	}
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
}

function setDirectory()
{
	dir = getDirectory("Choose a Source Directory");
	file_list = getFileList(dir);
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir))File.makeDirectory(output_dir);
	log_path = output_dir+"Log.txt";
}

function setFileNames(prefix)
{
	nuclei_roi_set 		= output_dir+prefix+"_nuclei_roi_set.zip";
	nuclei_results 		= output_dir+prefix+"_nuclei_results.txt";
	cells_roi_set 		= output_dir+prefix+"_cells_roi_set.zip";
	cells_results		= output_dir+prefix+"_cells_results.txt";
	folds_roi_set 		= output_dir+prefix+"_folds_roi_set.zip";
	folds_results		= output_dir+prefix+"_folds_results.txt";
	spots_a_roi_set		= output_dir+prefix+"_spots_a_roi_set.zip";
	spots_b_roi_set		= output_dir+prefix+"_spots_b_roi_set.zip";
	spots_a_results		= output_dir+prefix+"_spots_a_results.txt";
	spots_b_results		= output_dir+prefix+"_spots_b_results.txt";
	plaques_roi_set 	= output_dir+prefix+"_plaques_roi_set.zip";
	plaques_results		= output_dir+prefix+"_plaques_results.txt";
	ind_plaques_results	= output_dir+prefix+"_individual_plaques_results.txt";
	distance_results	= output_dir+prefix+"_distance_results.txt";
	results				= output_dir+prefix+"_summary.txt";
}

function scanFiles()
{
	prefixes = newArray(0);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,suffix) && indexOf(path,"flatfield")<0)
		{
			print(path);
			prefixes = Array.concat(prefixes,substring(file_list[i],0,lastIndexOf(file_list[i],suffix)));			
		}
	}
	return prefixes;
}

function setup()
{
	setOptions();
	Dialog.createNonBlocking("CellBlocks_v22 Settings");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("------------------------------------   General parameters  ------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Define here which regions you wish to analyze and whether you want to perform a distance calculation between ROI sets\n(not all comparisons are currently available)\n", 12, "#999999");
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Image Type", file_types, suffix);
	Dialog.addNumber("Pixel Size", pixel_size, 3, 5, micron+"");
	Dialog.addToSameRow();
	Dialog.addNumber("Field Number",fields, 0, 5, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Channel Number", channels, 0, 5, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Detectable Objects:");
	labels = newArray(5);defaults = newArray(5);
	labels[0] = "Nuclei";		defaults[0] = segment_nuclei;
	labels[1] = "Cells";		defaults[1] = segment_cells;
	labels[2] = "Spots";		defaults[2] = segment_spots;
	labels[3] = "Folds";		defaults[3] = segment_folds;
	labels[4] = "Plaques";		defaults[4] = segment_plaques;
	Dialog.setInsets(0,150,0);
	Dialog.addCheckboxGroup(1,5,labels,defaults);
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Measure Distances",measure_distances);
	Dialog.addToSameRow();
	Dialog.addChoice("From (secondary):", objects, secondary_objects);
	Dialog.addToSameRow();
	Dialog.addChoice("To/In (primary):", objects, primary_objects);
	Dialog.addToSameRow();
	Dialog.addCheckbox("In (off=to)",distance_within);
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Texture Analysis",texture_analysis);
	Dialog.setInsets(20,0,0);
	Dialog.addMessage("------------------------------------   Nuclei parameters  ------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Nuclei are segmented by classic smoothing or Laplacian enhancement, or using a trained classifier (Stardist).\nIf the latter is used, standard settings are not considered, but object filtering is always applied.\n", 12, "#999999");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Nuclei Channel", nuclei_channel, 0, 4, "");
	Dialog.addToSameRow();
	Dialog.addChoice("Segmentation Method", nuclei_segmentation_methods, nuclei_segmentation_method);
	Dialog.addToSameRow();
	Dialog.addNumber("Blur Radius", nuclei_filter_scale, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Standard Settings:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addCheckbox("Background Subtraction", nuclei_background);
	Dialog.addToSameRow();
	Dialog.addCheckbox("Contrast Enhancement", nuclei_clahe);
	Dialog.setInsets(0,0,0);	
	Dialog.addCheckbox("Watershed Separation", nuclei_watershed);
	Dialog.addToSameRow();
	Dialog.addChoice("Threshold Method", threshold_methods, nuclei_threshold);
	Dialog.addToSameRow();
	Dialog.addNumber("Fixed Threshold", nuclei_fixed_threshold_value, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Stardist settings:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Probability", nuclei_probability, 2, 4, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Tolerated overlap", nuclei_overlap, 2, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Object filters:\n");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Min. Circularity", nuclei_min_circularity, 2, 5, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Min. Area", nuclei_min_area, 0, 5, micron+"2");
	Dialog.addToSameRow();
	Dialog.addNumber("Max. Area", nuclei_max_area, 0, 5, micron+"2");	
	Dialog.setInsets(20,0,0);
	Dialog.addMessage("------------------------------------   Cell parameters  ------------------------------------", 14, "#ff0000");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Cells are segmented by thresholding a selected channel of a cell stain or if absent, by dilating nuclear ROI seeds.\nRegion growing is performed for a defined number of iterations, or if iterations is set to 0, sheer Voronoi tesselation is applied. \n", 12, "#999999");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Cell Channel", cells_channel,0,4," ");
	Dialog.addToSameRow();
	Dialog.addChoice("Segmentation Method", cells_segmentation_methods, cells_segmentation_method);
	Dialog.addToSameRow();
	Dialog.addCheckbox("Exclude Nuclei",exclude_nuclei);
	Dialog.addToSameRow();
	Dialog.addNumber("Filter scale ",cells_filter_scale,0,3,"");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------------------------------------------------------------------------------- ", 14, "#dddddd");
	Dialog.setInsets(0,0,0);	
	Dialog.addMessage("Threshold settings:                                                                                                      Dilation settings:                Cellpose settings:\n");

	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Threshold Method", threshold_methods, cells_threshold);
	Dialog.addToSameRow();
	Dialog.addNumber("Fixed Threshold", cells_fixed_threshold_value, 0, 4, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Grow cycles", iterations, 0, 5, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Expected diameter", cells_diameter, 0, 5, "");
	Dialog.show();

	print("%%%%%%%%%%%%%%%%CellBlocks_V22%%%%%%%%%%%%%%%%%%%%%%");
	suffix							= Dialog.getChoice();		print("Image Type:",suffix);
	pixel_size						= Dialog.getNumber(); 		print("Pixel Size:",pixel_size);
	fields 							= Dialog.getNumber();	 	print("Fields:",fields);
	channels 						= Dialog.getNumber();		print("Channels:",channels);
	segment_nuclei					= Dialog.getCheckbox(); 	print("Segment Nuclei:",segment_nuclei);
	segment_cells 					= Dialog.getCheckbox();		print("Segment Cells:",segment_cells);
	segment_spots 					= Dialog.getCheckbox();		print("Segment Spots:",segment_spots);
	segment_folds					= Dialog.getCheckbox();		print("Segment Folds:",segment_folds);
	segment_plaques					= Dialog.getCheckbox();		print("Segment Plaques:",segment_plaques);
	measure_distances				= Dialog.getCheckbox();		print("Measure Distances:",measure_distances);
	secondary_objects				= Dialog.getChoice();		print("Secondary Objects:",secondary_objects);
	primary_objects					= Dialog.getChoice();		print("Primary Objects:",primary_objects);
	distance_within					= Dialog.getCheckbox();		print("Measure within:",distance_within);	
	texture_analysis				= Dialog.getCheckbox();		print("Texture Analysis:",texture_analysis);
	nuclei_channel 					= Dialog.getNumber();		print("Nuclear Channel:",nuclei_channel);
	nuclei_segmentation_method		= Dialog.getChoice();		print("Nuclei Segmentation Method:",nuclei_segmentation_method);
	nuclei_filter_scale				= Dialog.getNumber();		print("Nuclei Filter Scale:",nuclei_filter_scale);
	nuclei_background				= Dialog.getCheckbox();		print("Background Subtraction:",nuclei_background);
	nuclei_clahe					= Dialog.getCheckbox();		print("Clahe:",nuclei_clahe);
	nuclei_watershed 				= Dialog.getCheckbox();		print("Watershed:",nuclei_watershed);
	nuclei_threshold				= Dialog.getChoice();		print("Nuclear Autothreshold:",nuclei_threshold);
	nuclei_fixed_threshold_value	= Dialog.getNumber();		print("Fixed Threshold Value:",nuclei_fixed_threshold_value);
	nuclei_probability 				= Dialog.getNumber();		print("Probability:",nuclei_probability);
	nuclei_overlap 					= Dialog.getNumber();		print("Overlap:",nuclei_overlap);
	nuclei_min_circularity			= Dialog.getNumber();		print("Min Circ:",nuclei_min_circularity);
	nuclei_min_area					= Dialog.getNumber();		print("Min Nuclear Size:",nuclei_min_area);
	nuclei_max_area					= Dialog.getNumber();		print("Max Nuclear Size:",nuclei_max_area);
	cells_channel 					= Dialog.getNumber();		print("Cell Channel:",cells_channel);
	cells_segmentation_method		= Dialog.getChoice();		print("Cell Segmentation Method:",cells_segmentation_method);
	if(cells_segmentation_method=="Trained Model")
	{
		model_dir = getDirectory("Where is the trained cell model?");
		list  = getFileList(model_dir);
		for(i=0;i<list.length;i++)
		{
			path = model_dir+list[i];
			if (endsWith(list[i], ".model"))
			{
				model = path;
				print("location of the cell classification model:",model);
			}			
		}
	}
	exclude_nuclei					= Dialog.getCheckbox();		print("Exclude Nuclear Area From Cell Analysis", exclude_nuclei);
	cells_filter_scale				= Dialog.getNumber();		print("Cell Filter Scale:", cells_filter_scale);
	cells_threshold					= Dialog.getChoice();		print("Cell Autothreshold:",cells_threshold);
	cells_fixed_threshold_value		= Dialog.getNumber();		print("Fixed Threshold Value:",cells_fixed_threshold_value);
	iterations						= Dialog.getNumber();		print("Region Growing Iterations:",iterations);
	cells_diameter					= Dialog.getNumber();		print("Cellpose diameter:",cells_diameter);
	
	
	if(segment_folds)
	{
		Dialog.createNonBlocking("CellBlocks Settings");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------   Folds parameters  ------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Detection of nuclear folds in lamin staining done by Laplacian enhancement and nuclear border erosion\n", 12, "#999999");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Folds Channel", folds_channel, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Segmentation Method", folds_segmentation_methods, folds_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Filter Scale",folds_filter_scale);
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Autothreshold Method",threshold_methods,folds_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold",folds_fixed_threshold_value);
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Border Shrink Size",shrink);
		Dialog.addToSameRow();
		Dialog.addNumber("Min. Fold Size",folds_min_area);	
		Dialog.setInsets(0,0,0);	
		Dialog.addCheckbox("Skeletonize Folds",skeletonize_folds);
		Dialog.show();

		folds_channel				= Dialog.getNumber();		print("Folds Channel:",folds_channel);
		folds_segmentation_method	= Dialog.getChoice();		print("Segmentation Method:",folds_segmentation_method);
		folds_filter_scale 			= Dialog.getNumber();		print("Filter scale:",folds_filter_scale);
		folds_threshold 			= Dialog.getChoice();		print("Spot AutoThreshold:",folds_threshold);
		folds_fixed_threshold_value	= Dialog.getNumber();		print("Fixed Threshold:",folds_fixed_threshold_value);
		shrink			 			= Dialog.getNumber();		print("Shrink Size:", shrink);
		folds_min_area 				= Dialog.getNumber();		print("Min. Fold Size:",folds_min_area);
		skeletonize_folds 			= Dialog.getCheckbox();		print("Skeletonize Folds:", skeletonize_folds);

		print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	}
	
	if(segment_spots)
	{
		Dialog.createNonBlocking("CellBlocks_v22 Settings");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------   Spot parameters  ------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Detection of spot-like structures in 1 or 2 channels by thresholding after Laplace/Gauss enhancement\nIf there is only one spot channel, set the second to 0. If both are present, colocalize allows detecting reciprocal overlap.\nWhen applying a fixed threshold in combination with laplace/multi-scale enhancement, negative threshold values should be used. ", 12, "#999999");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Spot Channel A", spots_a_channel, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Segmentation Method", spot_segmentation_methods, spots_a_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Filter Scale", spots_a_filter_scale, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Method", threshold_methods, spots_a_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_a_fixed_threshold_value, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Spot Size", spots_a_min_area, 0, 4, "pixels");
		Dialog.addToSameRow();
		Dialog.addNumber("Max. Spot Size", spots_a_max_area, 0, 4, "pixels");
		Dialog.setInsets(20,0,0);
		Dialog.addNumber("Spot Channel B", spots_b_channel, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Segmentation Method", spot_segmentation_methods, spots_b_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Filter Scale", spots_b_filter_scale, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Method", threshold_methods, spots_b_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_b_fixed_threshold_value, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Spot Size", spots_b_min_area, 0, 4, "px");
		Dialog.addToSameRow();
		Dialog.addNumber("Max. Spot Size", spots_b_max_area, 0, 4, "px");
		Dialog.setInsets(20,0,0);
		Dialog.addCheckbox("Colocalize Spots",colocalize_spots);
		Dialog.show();
		
		spots_a_channel					= Dialog.getNumber();		print("Spot Channel A:",spots_a_channel);
		spots_a_segmentation_method		= Dialog.getChoice();		print("Spot Segmentation Method:",spots_a_segmentation_method);
		spots_a_filter_scale 			= Dialog.getNumber();		print("Spot Filter Size:",spots_a_filter_scale);
		spots_a_threshold 				= Dialog.getChoice();		print("Spot AutoThreshold:",spots_a_threshold);
		spots_a_fixed_threshold_value  	= Dialog.getNumber();		print("Fixed Threshold Value:",spots_a_fixed_threshold_value);
		spots_a_min_area	 			= Dialog.getNumber();		print("Min. spot size:", spots_a_min_area);
		spots_a_max_area 				= Dialog.getNumber();		print("Max. spot size:",spots_a_max_area);
		
		spots_b_channel					= Dialog.getNumber();		print("Spot Channel B:",spots_b_channel);
		spots_b_segmentation_method		= Dialog.getChoice();		print("Spot Segmentation Method:",spots_b_segmentation_method);
		spots_b_filter_scale 			= Dialog.getNumber();		print("Spot Filter Scale:",spots_b_filter_scale);
		spots_b_threshold 				= Dialog.getChoice();		print("Spot AutoThreshold:",spots_b_threshold);
		spots_b_fixed_threshold_value  	= Dialog.getNumber();		print("Fixed Threshold Value:",spots_b_fixed_threshold_value);
		spots_b_min_area	 			= Dialog.getNumber();		print("Min. Spot Size:", spots_b_min_area);
		spots_b_max_area 				= Dialog.getNumber();		print("Max. Spot Size:",spots_b_max_area);
		
		colocalize_spots 				= Dialog.getCheckbox();		print("Colocalize spot channels",colocalize_spots); 
		
		print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	}

	if(segment_plaques)
	{
		Dialog.createNonBlocking("CellBlocks_v22 Settings");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------   Plaque parameters  ------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Detection of plaques that contain multiple nuclei by gross image thresholding (w or wo preprocessing).\n", 12, "#999999");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Plaque Channel", plaques_channel, 0, 4, "");
		Dialog.setInsets(10,0,0);
		Dialog.addCheckbox("Background Subtraction", plaques_background);
		Dialog.addToSameRow();
		Dialog.addCheckbox("Contrast Enhancement", plaques_clahe);
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Gaussian Scale", plaques_filter_scale, 0, 4, "");
		Dialog.addToSameRow();
		Dialog.addChoice("Threshold Method", threshold_methods, plaques_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", plaques_fixed_threshold_value, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Plaque Size", plaques_min_area, 0, 4, "px");
		Dialog.addToSameRow();
		Dialog.addNumber("Max. Plaque Size", plaques_max_area, 0, 4, "px");
		Dialog.setInsets(0,0,0);
		Dialog.addCheckbox("Watershed", plaques_watershed);
		Dialog.show();
		
		plaques_channel					= Dialog.getNumber();		print("Plaque Channel:",plaques_channel);
		plaques_background				= Dialog.getCheckbox();		print("Background Subtraction:", plaques_background);
		plaques_clahe					= Dialog.getCheckbox();		print("CLAHE:", plaques_clahe);
		plaques_filter_scale 			= Dialog.getNumber();		print("Plaque Filter Size:",plaques_filter_scale);
		plaques_threshold 				= Dialog.getChoice();		print("AutoThreshold Method:",plaques_threshold);
		plaques_fixed_threshold_value  	= Dialog.getNumber();		print("Fixed Threshold Value:",plaques_fixed_threshold_value);
		plaques_min_area	 			= Dialog.getNumber();		print("Min. Plaque Size:", plaques_min_area);
		plaques_max_area 				= Dialog.getNumber();		print("Max. Plaque Size:",plaques_max_area);
		plaques_watershed				= Dialog.getCheckbox();		print("Watershed:", plaques_watershed);
				
		print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	}
}

function calibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!=micron && unit!="microns" && unit!="micron")run("Properties...", " unit="+micron+" pixel_width="+pixel_size+" pixel_height="+pixel_size);
	else pixel_size = pixelWidth;
}

function decalibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!="pixel")run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
}

function segmentNuclei(id,c,sel)
{
	// input = multichannel image, output = roiset of nuclear ROIs and if(sel==1) mask incl. border objects
	// output = an image (mid) that contains all ROIs (also touching borders) and roiset of full nuclei
	mid = 0;
	selectImage(id);
	image_width = getWidth;
	image_height = getHeight;
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID; // the nuclear channel image to be turned into a binary image
	calibrateImage(cid);
	// preprocess the image
	selectImage(cid);
	if(nuclei_clahe)run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	if(nuclei_background)run("Subtract Background...", "rolling="+round(30/pixel_size));
	if(nuclei_segmentation_method != "Laplace" && nuclei_filter_scale > 0)run("Gaussian Blur...", "sigma="+nuclei_filter_scale);
	else if(nuclei_segmentation_method == "Laplace")
	{
		run("FeatureJ Laplacian", "compute smoothing="+nuclei_filter_scale); // scale to be adapted depending on nuclei size and SNR
		selectImage(cid); close;
		selectImage("copy Laplacian");
		rename("copy");
		cid = getImageID;
		selectImage(cid);
	}
	if(nuclei_segmentation_method != "Stardist")
	{
		if(nuclei_threshold=="Fixed")
		{
			if(nuclei_segmentation_method == "Laplace")
			{
				setAutoThreshold("Default ");
				getThreshold(mit,mat); 
				setThreshold(mit,nuclei_fixed_threshold_value);
			}
			else 
			{
				setAutoThreshold("Default dark");
				getThreshold(mit,mat); 
				setThreshold(nuclei_fixed_threshold_value, mat);
			}
		}
		else {
			if(nuclei_segmentation_method == "Laplace")setAutoThreshold(nuclei_threshold); 
			else setAutoThreshold(nuclei_threshold+" dark"); 
		}
		getThreshold(mit,mat); print("Nuclei Threshold:",mit,mat);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Fill Holes");
		if(nuclei_watershed)run("Watershed");
	}
	else if(nuclei_segmentation_method == "Stardist")
	{
		selectImage(cid);
		run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
		+"args=['input':"+'copy'+", 'modelChoice':'Versatile (fluorescent nuclei)',"
		+"'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99',"
		+"'probThresh':'"+nuclei_probability+"', 'nmsThresh':'"+nuclei_overlap+"', 'outputType':'ROI Manager', 'nTiles':'1', "
		+"'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', "
		+"'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
		selectImage(cid); close;
		newImage("copy", "8-bit black", image_width, image_height, 1);
		calibrateImage(cid);
		cid = getImageID;
		selectImage(cid);
		selectImage("copy");
		nr_rois = roiManager("count");
		for(r = 0; r < nr_rois; r++)
		{
			roiManager("select",r);
			run("Enlarge...", "enlarge=1");
			run("Clear");
			run("Enlarge...", "enlarge=-1");
			run("Fill");
		}
		roiManager("Deselect");
		roiManager("reset");
		setThreshold(1,255);
		run("Convert to Mask");
	}
	run("Set Measurements...", "area centroid mean integrated redirect=None decimal=2");
	if(sel)
	{
		selectImage(cid);
		run("Analyze Particles...", "size="+nuclei_min_area+"-"+nuclei_max_area+" circularity="+nuclei_min_circularity+"-1.00 show=Masks clear include");
		if(isOpen("Mask of copy"))
		{
			selectWindow("Mask of copy"); 
			mid = getImageID;   // the full mask of all particles (for more accurate cell segmentation)
		}
	}
	selectImage(cid);
	// only apply rim to the nuclei you actually want to analyze
	if(rim>0)run("Specify...", "width="+(image_width-image_width*rim/100)+" height="+(image_height-image_height*rim/100)+" x="+image_width*rim/200+" y="+image_height*rim/200);
	run("Analyze Particles...", "size="+nuclei_min_area+"-"+nuclei_max_area+" circularity="+nuclei_min_circularity+"-1.00 show=Nothing exclude clear include add");
	rmc = roiManager("count"); print(rmc,"Nuc. ROI");
	if(rmc==0 && isOpen(mid)){selectImage(mid); close; mid=0;}
	for(i=0;i<rmc;i++)
	{
		roiManager("select",i);
		if(i<9)roiManager("Rename","000"+i+1);
		else if(i<99)roiManager("Rename","00"+i+1);	
		else if(i<999)roiManager("Rename","0"+i+1);	
		else roiManager("Rename",i+1);
	}
	run("Clear Results");
	roiManager("deselect"); 
	roiManager("Measure");
	selectImage(cid); close; 
	return mid;
}

function flatfield_correct(id){
	//	flatField does a flatfield field correction, using an image of fluorescent plastic
	//	if flatfield contains different channels and input image only one, c directs to the correct channel to select in the flatfield stack
	selectImage(id);
	title = getTitle;
	flatpath = dir+"flatfield.tif";
	open(flatpath);
	fid = getImageID;
	selectImage(fid);
	getRawStatistics(n,flatmean);
	imageCalculator("Divide create 32-bit stack", id,fid); 
	selectImage(fid); close;
	selectImage(id); close;
	selectWindow("Result of "+title);
	id = getImageID;
	selectImage(id);
	rename(title);  
	run("Multiply...","value="+flatmean+" stack");
	return id;
}

function segmentRegions(id, mid, c, iterations){
	selectImage(mid); 
	run("Select None");
	name = getTitle;
	
	// generate voronoi regions from all detected nuclei (including edges) to have some rough boundaries between touching cells
	run("Duplicate...","title=voronoi");
	vid = getImageID;
	selectImage(vid);
	run("Voronoi");
	setThreshold(1, 255);
	run("Convert to Mask");
	run("Invert");
	
	// generate dilated nuclei (using more accurate EDM mapping (by x iterations) requires the biovioxxel package) - now just using the enlarge function
	if(cells_segmentation_method=="Dilation" && iterations != 0)
	{
		selectImage(mid);
		run("Duplicate...","title=dilate");
		did = getImageID;
		selectImage(did); 
		//run("EDM Binary Operations", "iterations="+iterations+" operation=dilate"); 
		run("Create Selection"); 
		run("Enlarge...", "enlarge="+iterations+" pixel");
		roiManager("Add");
		run("Select All");
		run("Fill");
		sel = roiManager("count")-1;
		roiManager("select", sel);
		run("Clear", "slice");
		roiManager("Delete");
		run("Select None");
		run("Invert LUT");
		imageCalculator("AND create", "voronoi","dilate");
		selectImage(did); close; 
		selectImage(vid); close;
		selectWindow("Result of voronoi");
		if(!is("Inverting LUT"))run("Invert LUT");
		//run("Invert");
		vid = getImageID;
		selectImage(vid);
		rename("Cell_ROI");
	}
	
	// voronoi
	if(cells_segmentation_method=="Dilation" && iterations == 0)
	{
		selectImage(vid);
		rename("Cell_ROI");
	}

	// use a cellular counterstain to make accurate cell ROIs
	if(cells_segmentation_method!="Dilation"){
		selectImage(id);
		run("Select None");
		if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
		else{setSlice(c);run("Duplicate...","title=copy ");}
		cid = getImageID;
		selectImage(cid);
		if(cells_filter_scale > 0)run("Gaussian Blur...", "sigma="+cells_filter_scale);
		if(cells_segmentation_method=="Trained Model")
		{
			print("Applying trained model, takes time!"); //model needs to be in the FIJI/scripts foldr along with the Trained_Glia_Segmentation.bsh script in teh FIJI/scripts folder
			run("Trained Cell Segmentation","model="+model+" image="+cid);
			selectImage(cid); close;
			selectWindow("Classification result");
			rename("copy");
			cid = getImageID();		
			setAutoThreshold("Default "); 	
		}
		else if(cells_segmentation_method=="Threshold")
		{
			if(cells_threshold=="Fixed"){
				setAutoThreshold("Default dark");
				getThreshold(mit,mat); 
				setThreshold(cells_fixed_threshold_value,mat);
			}
			else {
				setAutoThreshold(cells_threshold+" dark"); 
			}
		}
		else if(cells_segmentation_method=="Cellpose"){
			print("Cellpose Start");

			//Create RGB image for cellpose: R=1, G=2, B=3
			selectImage(id);
			run("Select None");
			if(Stack.isHyperstack)run("Duplicate...", "title=Nuclei_Red duplicate channels="+nuclei_channel);	
			else{ setSlice(nuclei_channel); run("Duplicate...","title=Nuclei_Red ");}
			redid=getImageID();
			selectImage(id);
			run("Select None");
			if(Stack.isHyperstack)run("Duplicate...", "title=Cells_Green duplicate channels="+cells_channel);	
			else{setSlice(cells_channel);run("Duplicate...","title=Cells_Green ");}
			greenid=getImageID();
			run("Merge Channels...", "c1=[Nuclei_Red] c2=[Cells_Green] create");
			compid=getImageID();
			run("RGB Color");
			rgbid=getImageID();
			selectImage(cid); close;

			//Run Cellpose and exclude roi borders from label map
			run("Cellpose Script Wrapper", "imp="+id+" diameter="+cells_diameter+" channel_nuc="+1+" channel_cyto="+2);
			selectWindow("Cellpose_Label");
			cid = getImageID();
			rename("copy");
			selectImage(rgbid); close;
			selectImage(compid); close;
			setForegroundColor(0, 0, 0);
			getRawStatistics(nPixels, mean, min, max, std, histogram);	
			maxMask=max;
			minPixels=(nuclei_min_area/(pixel_size*pixel_size))/10;  //min pixel number equal to 1/10 of min nuc area	
			maxPixels=nPixels;		
			for(i=1; i<=maxMask; i++)
			{
				selectImage(cid);
				setThreshold(i, i, "raw");
				run("Create Selection");
				getRawStatistics(nPixels, mean, min, max, std, histogram);	
				
				if(nPixels > minPixels && nPixels!=maxPixels){
					run("Draw", "slice");
				}
			}
			setThreshold(1, maxMask, "raw");
			setForegroundColor(255, 255, 255);
			print("Cellpose End");
		}
		selectImage(cid);
		getThreshold(mit,mat); 
		print("Cell Threshold:",mit,mat);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		//run("Fill Holes");	
		
		if(cells_segmentation_method!="Cellpose")
		{
			imageCalculator("AND create", "voronoi","copy"); // apply the voronoi bounderies to the cell ROIs
			selectImage(cid); close; 
			selectImage(vid); close;
			selectWindow("Result of voronoi");
			vid = getImageID;
		}else { //cellpose
			selectImage(vid); close;
			selectImage(cid);
			vid=getImageID();
		}
		selectImage(vid);
		rename("Cell_ROI");
	}
	
	// keep only the non-excluded nuclei (pos nuclei)
	newImage("posnuclei", "16-bit Black", image_width, image_height, 1); 
	pid = getImageID;
	selectImage(pid); 
	roiManager("Deselect");
	roiManager("Fill");
	run("Convert to Mask");
	
	// fuse cell and nuclei image to avoid non-overlapping ROI (not with cellpose)
	if(cells_segmentation_method!="Cellpose"){
		imageCalculator("OR create", "Cell_ROI","posnuclei"); 
		selectImage(vid); close;
		selectWindow("Result of Cell_ROI");
		rename("Cells");
		vid = getImageID;
	}else{
		selectImage(vid);
		rename("Cells");
	}
	
	// make labeled nuclei centroid mask
	newImage("centroid", "16-bit Black", image_width, image_height, 1); 
	pidc=getImageID();
	calibrateImage(pidc);
	run("16-bit"); //Needed if > 255 nuclei detected
	rmc = roiManager("count"); // rmc = number of nuclei
	print(rmc,"retained nuclei");
	selectImage(pidc);
	run("Set Measurements...", "centroid redirect=None decimal=4");
	roiManager("measure");
	for(i=0;i<rmc;i++)
	{
		Mx=getResult("X", i)/pixel_size;
		My=getResult("Y", i)/pixel_size;
	 	makeOval(Mx-2, My-2, 4, 4);
		run("Set...", "value="+i+1);		
	}
	run("Select None");
	
	// Add cell rois and rename with matching nuclear index 
	selectImage(vid);
	run("Analyze Particles...", "size="+nuclei_min_area+"-Infinity circularity=0.00-1.00 show=Nothing add");
	rmcb = roiManager("count"); // number of cell regions larger than a nucleus 
	print(rmcb-rmc,"Number of detected cell regions larger than a min. nuclear area"); 
	selectImage(pidc);
	for(i=rmcb-1;i>=rmc;i--)
	{
		roiManager("select",i);
		getRawStatistics(np,mean,min,max);
		if(max==0){roiManager("delete");} 				// no nuc so not retained
		else if(max<10)roiManager("Rename","000"+max);	// assigned to correct nuc
		else if(max<100)roiManager("Rename","00"+max);	
		else if(max<1000)roiManager("Rename","0"+max);	
		else roiManager("Rename",max);
	}	
	selectImage(pidc); close();
	
	rmcc = roiManager("count"); //all cell regions with one nucleus	
	print(rmcc-rmc,"Number of unique cell regions that overlap with a nucleus"); 
	
	// exclude nuclei
	if(rmcc>rmc)
	{
		if(exclude_nuclei) //define cytoplasmic regions (cells without nuclei)
		{
			//Each Nucleus has a cell region, since the nuclei and cell image were merged
			//Check for cells that are equal in size to nucleus --> No cytoplasm --> Do not calculate XOR --> No Cytoplasm region added --> 0 in summary file
			roiManager("Sort");
			index = 0;
			for(i=0;i<rmc;i++)
			{
				index=i*2;
				roiManager("select",index); 
				roi_name_a = Roi.getName();
				getRawStatistics(np_a); //area
				roiManager("select", index+1); 
				roi_name_b = Roi.getName(); 
				getRawStatistics(np_b); //area
				if (matches(roi_name_a, roi_name_b) && np_a!=np_b)
				{ 
					couple = newArray(index,index+1);
					roiManager("select",couple); 
					roiManager("XOR"); 
					roiManager("Add");
					roiManager("select",roiManager("count")-1);
					roiManager("Rename",roi_name_a);
				}else if(matches(roi_name_a, roi_name_b) && np_a==np_b)
				{
					print("Matched with Cell ROI",roi_name_b,"But discarded due to equal size (no cytoplasm detected)");
				}
			}
			roiManager("select",Array.getSequence(rmcc));
			roiManager("Delete"); 			
		} 
		else if(!exclude_nuclei)
		{
			roiManager("select",Array.getSequence(rmc)); 
			roiManager("Delete"); //remove all nuclei ROIs and retain cell ROIs
		} 
		run("Select None");
		cell_nr = roiManager("count");		
	}
	else 
	{
		cell_nr = 0;
	}
	print(cell_nr, "Number of unique cell regions retained with area larger than a nucleus");
	selectImage(mid); close;
	selectImage(vid); close;
	selectImage(pid); close;
	if(isOpen("copy")){selectWindow("copy");close;}
	return cell_nr;
}

function segmentFolds(id,c)
{
	//	segments the prominent lamina structures (inside and peripheral)
	selectImage(id);
	run("Select None");
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID;
	selectImage(cid);
	decalibrateImage(cid);
	title = getTitle;
	if(folds_segmentation_method=="Gauss")
	{
		run("Duplicate...","title=folds");
		run("Gaussian Blur...", "sigma="+folds_filter_scale);
		lid = getImageID;
	}
	else if(folds_segmentation_method=="Sobel")
	{
		run("FeatureJ Edges", "compute smoothing="+folds_filter_scale);
		lid = getImageID;
		selectImage(lid); rename("folds");
	}
	else if(folds_segmentation_method=="Laplace")
	{
		run("FeatureJ Laplacian", "compute smoothing="+folds_filter_scale);
		lid = getImageID;
		selectImage(lid); rename("folds");
	}
	selectImage(cid); close;
	selectImage(lid);
	if(folds_threshold=="Fixed")
	{
		if(folds_segmentation_method!="Laplace")setAutoThreshold("Default dark ");
		else setAutoThreshold("Default ");
		getThreshold(mit,mat); 
		print("Original threshold settings:",mit,mat);
		if(folds_segmentation_method!="Laplace")setThreshold(folds_fixed_threshold_value,maxOf(mat,folds_fixed_threshold_value));
		else setThreshold(minOf(mit,-folds_fixed_threshold_value),-folds_fixed_threshold_value);
	}
	else 
	{
		if(folds_segmentation_method!="Laplace")setAutoThreshold(folds_threshold+" dark ");
		else setAutoThreshold(folds_threshold+" ");
	}
	getThreshold(mit,mat); 
	print("Threshold:",mit,mat);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	if(skeletonize_folds)run("Skeletonize");
	run("Analyze Particles...", "size=0-Infinity pixel circularity=0-1.00 show=Masks");
	selectImage(lid); close;
	selectWindow("Mask of folds");
	fid = getImageID;
	//	trims ROIs down to what is supposed to the inner side (defined by enlarge)
	selectImage(fid);
	roiManager("Deselect");
	rmc = roiManager("count"); 
	if(rmc==1)roiManager("select",0);
	else roiManager("Combine");
	run("Enlarge...", "enlarge=-"+shrink+" pixel");
	run("Make Inverse");
	run("Fill");
	run("Select None");
	roiManager("reset");
	run("Analyze Particles...", "size="+folds_min_area+"-Infinity pixel circularity=0-1.00 show=Nothing add");
	fnr = roiManager("count");
	selectImage(fid); close;
	return fnr;
}

function segmentSpots(id,c,args)
{
	spot_channel				= args[0];
	spot_segmentation_method	= args[1];
	spot_threshold_method		= args[2];
	spot_fixed_threshold_value	= args[3];
	scale						= args[4];
	spot_min_area				= args[5];
	spot_max_area				= args[6];
	selectImage(id);
	run("Select None");
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID;
	decalibrateImage(cid);
	selectImage(cid);
	title = getTitle;
	//run("Enhance Local Contrast (nuclei_clahe)", "blocksize=125 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	if(spot_segmentation_method=="Gauss")
	{
		run("Duplicate...","title=Gauss");
		run("Gaussian Blur...", "sigma="+scale);
		lap = getImageID;
	}
	else if(spot_segmentation_method=="Laplace")
	{
		run("FeatureJ Laplacian", "compute smoothing="+scale);
		lap = getImageID;
	}
	else if(spot_segmentation_method=="Multi-Scale") 
	{
		e = 0;
		while(e<scale)
		{			
			e++;
			selectImage(cid);
			run("FeatureJ Laplacian", "compute smoothing="+e);
			selectWindow(title+" Laplacian");
			run("Multiply...","value="+e*e);
			rename("scale "+e);
			eid = getImageID;
			if(e>1)
			{
				selectImage(eid);run("Select All");run("Copy");close;
				selectImage(fid);run("Add Slice");run("Paste");
			}
			else fid = getImageID;
		}
		selectImage(fid);
		run("Z Project...", "start=1 projection=[Sum Slices]");
		lap = getImageID;
		selectImage(fid); close;
	}
	selectImage(lap);	
	if(spot_threshold_method=="Fixed")
	{
		if(spot_segmentation_method == "Laplace" || spot_segmentation_method == "Multi-Scale")
		{
			setAutoThreshold("Default ");
			getThreshold(mit,mat); 
			setThreshold(mit,spot_fixed_threshold_value);
		}
		else 
		{
			setAutoThreshold("Default dark");
			getThreshold(mit,mat); 
			setThreshold(spot_fixed_threshold_value, mat);
		}	
	}
	else 
	{
		if(spot_segmentation_method=="Gauss")setAutoThreshold(spot_threshold_method + " dark");
		else setAutoThreshold(spot_threshold_method+" ");
	}
	getThreshold(mit,mat); 
	print("Threshold:",mit,mat);
	run("Set Measurements...", "  area min mean redirect=["+title+"] decimal=4");
	run("Analyze Particles...", "size="+spot_min_area+"-"+spot_max_area+" circularity=0.00-1.00 show=Nothing display clear include add");
	snr = roiManager("count"); print(snr,"spots");
	if(snr>10000){print("excessive number, hence reset"); snr=0; roiManager("reset");}
	//to avoid excessive spot finding when there are no true spots
	selectImage(lap); close;
	selectImage(cid); close;
	return snr;
}

function segmentPlaques(id, c)
{
	// input = multichannel image, output = 1 ROI of segmented material
	sel = 0;
	selectImage(id); 
	run("Select None");
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID; 
	calibrateImage(cid);
	// preprocess the image
	selectImage(cid);
	if(plaques_clahe)run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	if(plaques_background)run("Subtract Background...", "rolling=100");
	if(plaques_filter_scale > 0)run("Gaussian Blur...", "sigma="+plaques_filter_scale);
	if(plaques_threshold=="Fixed")
	{
		setAutoThreshold("Default dark");
		getThreshold(mit,mat); 
		print("Original threshold settings:",mit,mat);
		setThreshold(plaques_fixed_threshold_value,mat);
	}
	else setAutoThreshold(plaques_threshold+" dark");
	getThreshold(mit,mat); 
	print("Threshold settings:",mit,mat);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	if(plaques_watershed)run("Watershed");
	run("Analyze Particles...", "size="+plaques_min_area+"-"+plaques_max_area+" show=Masks clear");
	selectWindow("Mask of copy"); 
	mid = getImageID;
	selectImage(mid);
	run("Invert"); 
	run("Create Selection");
	if(selectionType>-1)
	{
		roiManager("add");
		sel = 1;
	}
	else print("No plaques detected");
	selectImage(cid); close;
	selectImage(mid); close;
	return(sel);
}

function analyzeRegions(id)
{
	erase(0); 
	mask = 0;
	readout = 1;
	//	analyze cell rois
	selectImage(id);
	calibrateImage(id);
	if(File.exists(cells_roi_set))
	{
		run("Set Measurements...", "area mean standard modal min centroid center perimeter shape integrated median skewness kurtosis redirect=None decimal=4");
		roiManager("Open",cells_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}

		if(texture_analysis){
			textureAnalysis(id);
		}

		sortResults(); // organize results per channel
		
		if(!isOpen(mask))
		{
			newImage("Mask", "32-bit Black",image_width, image_height, 1); 	//	reference image for spot assignments
			mask = getImageID; 
		}
		selectImage(mask);
		for(j=0;j<rmc;j++)
		{
			roiManager("select",j);
			index = getInfo("roi.name");
			run("Set...", "value="+0-index);							//	negative values for cytoplasm, positive for nuclei
			setResult("Nuclei",j,index);
		}	
		saveAs("Measurements",cells_results);
		erase(0);
	}	
	//	analyze nuclear rois
	if(File.exists(nuclei_roi_set))
	{
		run("Set Measurements...", "area mean standard modal min centroid center perimeter shape feret's integrated median skewness kurtosis redirect=None decimal=4");
		roiManager("Open",nuclei_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		
		if(texture_analysis){
			textureAnalysis(id);
		}
	
		sortResults(); //organise results per channel
				
		if(!isOpen(mask))
		{
			newImage("Mask", "32-bit Black",image_width, image_height, 1); //	reference image for spot assignments 
			mask = getImageID;
		}
		selectImage(mask);
		for(j=0;j<rmc;j++)
		{
			roiManager("select",j);
			index = getInfo("roi.name");
			run("Set...", "value="+index);					//	negative values for cytoplasm, positive for nuclei
			setResult("Cell",j,index);
		}	
		run("Select None");
		updateResults;
		saveAs("Measurements",nuclei_results);
		erase(0);
	}	
		
	//	rudimentary colocalization analysis by binary overlap of spot ROIs requires bianry masks of spots
	if(colocalize_spots && File.exists(spots_a_roi_set) && File.exists(spots_b_roi_set))
	{
		roiManager("reset");
		roiManager("Open",spots_a_roi_set);
		selectImage(mask);
		run("Add Slice");
		setForegroundColor(255,255,255);
		setSlice(2);
		roiManager("Fill");
		roiManager("reset");
		roiManager("Open",spots_b_roi_set);
		selectImage(mask);
		run("Add Slice");
		setSlice(3);
		roiManager("Fill");
		run("Select None");
		roiManager("reset");
	}
	//	analyze spot rois
	if(File.exists(spots_a_roi_set))
	{
		selectImage(mask); 
		ms = nSlices;
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",spots_a_roi_set);
		spot_nr = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		// determine the location of the spots (cell vs. nucleus)
		selectImage(mask); setSlice(1);
		roiManager("deselect");
		roiManager("Measure");
		nindices = newArray(spot_nr);
		cindices = newArray(spot_nr);	
		
		for(j=0;j<spot_nr;j++)
		{
			min = getResult("Min",j);
			max = getResult("Max",j);
			if(max>0){
				nindices[j] = max; 
				if(exclude_nuclei){
					if(min > 0){cindices[j] = 0;}
					else{cindices[j] = -min;}
				}
				else{cindices[j] = max;}
			}else if(min<0){nindices[j] = 0; cindices[j] = -min;}
		}	
		run("Clear Results");
		// determine the colocalizing spots (one pixel overlap is sufficient)
		if(colocalize_spots && ms==3)
		{
			selectImage(mask); setSlice(3);
			roiManager("Measure");
			overlaps = newArray(spot_nr);
			for(j=0;j<spot_nr;j++)
			{
				max = getResult("Max",j);
				if(max>0){overlaps[j]=1;}
			}	
			selectWindow("Results"); run("Close");
		}
		IJ.renameResults("Temp","Results");
		for(j=0;j<spot_nr;j++)
		{
			if(colocalize_spots && ms==3)setResult("Coloc",j,overlaps[j]);
			setResult("Nucleus",j,nindices[j]);
			setResult("Cell",j,cindices[j]);
		}
		updateResults;
		saveAs("Measurements",spots_a_results);
		erase(0);
	}
	if(File.exists(spots_b_roi_set))
	{
		selectImage(mask); 
		ms = nSlices;
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",spots_b_roi_set);
		spot_nr = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		selectImage(mask);setSlice(1); 
		roiManager("deselect");
		roiManager("Measure");
		nindices = newArray(spot_nr);
		cindices = newArray(spot_nr);	
		for(j=0;j<spot_nr;j++)
		{
			min = getResult("Min",j);
			max = getResult("Max",j);
			if(max>0){
				nindices[j] = max; 
				if(exclude_nuclei){
					if(min > 0){cindices[j] = 0;}
					else{cindices[j] = -min;}
				}
				else{cindices[j] = max;}
			}else if(min<0){nindices[j] = 0; cindices[j] = -min;}
		}	
		run("Clear Results");
		// determine the colocalizing spots (one pixel overlap is sufficient)
		if(colocalize_spots && ms==3)
		{
			selectImage(mask); setSlice(2);
			roiManager("Measure");
			overlaps = newArray(spot_nr);
			for(j=0;j<spot_nr;j++)
			{
				max = getResult("Max",j);
				if(max>0){overlaps[j]=1;}
			}	
			selectWindow("Results"); run("Close");
		}
		IJ.renameResults("Temp","Results");
		for(j=0;j<spot_nr;j++)
		{
			if(colocalize_spots && ms==3)setResult("Coloc",j,overlaps[j]);
			setResult("Nucleus",j,nindices[j]);
			setResult("Cell",j,cindices[j]);
		}
		updateResults;
		saveAs("Measurements",spots_b_results);
		erase(0);
	}
	// fold rois
	if(File.exists(folds_roi_set))
	{
		selectImage(mask);
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",folds_roi_set);
		fold_nr = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		selectImage(mask); setSlice(1); 
		roiManager("deselect");
		roiManager("Measure");
		cindices = newArray(fold_nr);	
		nindices = newArray(fold_nr);
		for(j=0;j<fold_nr;j++)
		{
			min = getResult("Min",j);
			max = getResult("Max",j);
			if(max>0){
				nindices[j] = max; 
				if(exclude_nuclei){
					if(min > 0){cindices[j] = 0;}
					else{cindices[j] = -min;}
				}
				else{cindices[j] = max;}
			}else if(min<0){nindices[j] = 0; cindices[j] = -min;}
		}	
		run("Clear Results");
		IJ.renameResults("Temp","Results");
		for(j=0;j<fold_nr;j++)
		{
			setResult("Nucleus",j,nindices[j]);
			setResult("Cell",j,cindices[j]);
		}
		updateResults;
		saveAs("Measurements",folds_results);
		erase(0);

	}
	// analyze plaque rois
	if(segment_plaques)
	{
		selectImage(mask);
		run("Select None");
		getRawStatistics(image_area, mean, min, max);
		setThreshold (1,max);
		setOption("BlackBackground",false);
		run("Convert to Mask");
		getRawStatistics(image_area, image_mean);
		nuclear_area_image = image_area * image_mean/255;	
		nuclear_content_image = nuclear_area_image / image_area * 100;	
		run("Analyze Particles...","display");
		nuclei_number_image = nResults;	
		selectWindow("Results"); run("Close");
		if(File.exists(plaques_roi_set))
		{
			roiManager("Open",plaques_roi_set);
			selectImage(mask);
			roiManager("select",0);
			run("Analyze Particles...","display");
			nuclei_number_plaque = nResults;	
			selectWindow("Results"); run("Close");
			getRawStatistics(plaque_area, plaque_mean);	
		}
		else
		{
			plaque_area = 0;
			plaque_mean = 0;
			nuclei_number_plaque = 0;	
		}
		nuclei_number_ratio = nuclei_number_plaque / nuclei_number_image *100;	
		plaque_coverage = plaque_area / image_area * 100;
		nuclear_area_plaque = plaque_area * plaque_mean/255;		
		nuclear_content_plaque = nuclear_area_plaque / plaque_area * 100;		
		nuclear_content_ratio = nuclear_area_plaque / nuclear_area_image *100;
		setResult("Image_Area",0,image_area);
		setResult("Plaque_Area",0,plaque_area);
		setResult("Plaque_Coverage",0,plaque_coverage);
		setResult("Nuclei_Number_Image",0,nuclei_number_image);
		setResult("Nuclear_Area_Image",0,nuclear_area_image);
		setResult("Nuclear_Content_Image",0,nuclear_content_image);
		setResult("Nuclei_Number_Plaque",0,nuclei_number_plaque);
		setResult("Nuclear_Area_Plaque",0,nuclear_area_plaque);
		setResult("Nuclear_Content_Plaque",0,nuclear_content_plaque);
		setResult("Nuclear_Content_Ratio",0,nuclear_content_ratio);
		setResult("Nuclear_Number_Ratio",0,nuclei_number_ratio);
		updateResults();
		saveAs("Measurements",plaques_results);
		selectWindow("Results"); run("Close");
		if(File.exists(plaques_roi_set))
		{
			if(selectionType == 9) // only when composite
			{			
				roiManager("Select", 0);
				roiManager("Split");
				roiManager("Select", 0);
				roiManager("Delete");
			}
			selectImage(mask);
			roiManager("Deselect");
			run("Select None");
			run("Set Measurements...", "  area mean min redirect=None decimal=4");
			roiManager("Measure");
			for(r = 0; r < nResults; r++)
			{
				plaque_size = getResult("Area",r);
				nuclear_size = getResult("Mean",r)*plaque_size/255;
				plaque_nuclear_content = nuclear_size/plaque_size*100;
				setResult("Plaque_Area",r,plaque_size);
				setResult("Nuclear_Area",r,nuclear_size);
				setResult("Plaque_nuclear_content",r,plaque_nuclear_content);
			}
			updateResults();
			selectImage(mask); 
			run("Find Maxima...", "prominence=0 light output=[Single Points]");
			selectImage(mask); close;
			selectWindow("Mask Maxima"); 
			mask = getImageID;
			selectImage(mask);
			roi_count = roiManager("count");
			for(r = 0; r < roi_count; r++)
			{
				roiManager("select",r);
				getRawStatistics(plaque_area, nuclei_mean);	
				plaque_nuclei_nr = plaque_area*nuclei_mean / 255;
				setResult("count",r,plaque_nuclei_nr);
			}
			updateResults();
			saveAs("Measurements",ind_plaques_results);
			selectWindow("Results"); run("Close");		
		}
		erase(0);
	}
	if(isOpen(mask)){selectImage(mask); close;}
	else readout = 0;
	if(measure_distances)
	{
		run("Set Measurements...", "mean min redirect=None decimal=4");
		if(primary_objects=="Nuclei")primary_roi_set = nuclei_roi_set; 		
		else if(primary_objects=="Cells")primary_roi_set = cells_roi_set; 
		else if(primary_objects=="Folds")primary_roi_set = folds_roi_set; 
		else if(primary_objects=="Spots_a")primary_roi_set = spots_a_roi_set; 
		else if(primary_objects=="Spots_b")primary_roi_set = spots_b_roi_set; 	
		else primary_roi_set = plaques_roi_set; 
		if(secondary_objects=="Nuclei")secondary_roi_set = nuclei_roi_set; 		
		else if(secondary_objects=="Cells")secondary_roi_set = cells_roi_set; 
		else if(secondary_objects=="Folds")secondary_roi_set = folds_roi_set; 
		else if(secondary_objects=="Spots_a")secondary_roi_set = spots_a_roi_set; 
		else if(secondary_objects=="Spots_b")secondary_roi_set = spots_b_roi_set; 	
		else secondary_roi_set = plaques_roi_set; 
		if(File.exists(primary_roi_set) && File.exists(secondary_roi_set))
		{
			roiManager("Open",primary_roi_set);
			newImage("Mask", "8-bit Black",image_width, image_height, 1); 	
			mask = getImageID; 
			selectImage(mask);
			roiManager("Deselect");
			roiManager("Fill");
			roiManager("reset");
			run("Select None");
			if(distance_within)run("Invert");
			run("Options...", "iterations=1 count=1 edm=16-bit do=Nothing");
			run("Distance Map");
			selectImage(mask); 
			close;
			selectImage("EDM of Mask");
			mask = getImageID; 
			roiManager("Open",secondary_roi_set);
			selectImage(mask);
			roiManager("Measure");
			if(nResults>0)
			{
				//Table.renameColumn("Mean", "Mean.Distance");
				//Table.renameColumn("Max", "Max.Distance");
				//Table.renameColumn("Min", "Min.Distance");
				saveAs("Measurements",distance_results);
			}
			selectImage(mask); 
			close;
			erase(0);
		}
		else 
		{
			if(!File.exists(primary_roi_set))print("No distance calculation due to lack of primary ROI set");
			if(!File.exists(secondary_roi_set))print("No distance calculation due to lack of secondary ROI set");
		}	
	}
	return readout;
}

function textureAnalysis (id){
	IJ.renameResults("Results","Temp");
	run("Clear Results");
	rmc=roiManager("count");
	id=getImageID();
	run("Select None");
	
	angSecMom=newArray(0);
	invDiffMom=newArray(0);
	contrast=newArray(0);
	energy=newArray(0);
	entropy=newArray(0);
	homogeneity=newArray(0);
	variance=newArray(0);
	shade=newArray(0);
	prominence=newArray(0);
	inertia=newArray(0);
	correlation=newArray(0);
	
	var angles=newArray(0,45,90,135);
	for(c=1;c<=channels;c++){
		selectImage(id);
		if(Stack.isHyperstack){
			run("Duplicate...", "title=copy duplicate channels="+c);	
		}else{
			setSlice(c);
			run("Duplicate...","title=copy ");
		}

		resetMinAndMax;
		run("8-bit");
		idCh=getImageID();

		for(a=0;a<angles.length;a++){
			for(i=0; i<rmc;i++){
				selectImage(idCh);
				roiManager("select", i);
				run("GLCM Texture3", "enter=1 select="+angles[a]+" symmetrical angular contrast correlation inverse entropy energy inertia homogeneity prominence variance shade");
			}
			
		}
		selectImage(idCh);close;
	}
		
	for(i=0;i<nResults;i++){
		angSecMom=Array.concat(angSecMom,getResult("Angular Second Moment", i));
		invDiffMom=Array.concat(invDiffMom,getResult("Inverse Difference Moment", i));
		contrast=Array.concat(contrast,getResult("Contrast", i));
		energy=Array.concat(energy,getResult("Energy", i));
		entropy=Array.concat(entropy,getResult("Entropy", i));
		homogeneity=Array.concat(homogeneity,getResult("Homogeneity", i));
		variance=Array.concat(variance,getResult("Variance", i));
		shade=Array.concat(shade,getResult("Shade", i));
		prominence=Array.concat(prominence,getResult("Prominence", i));
		inertia=Array.concat(inertia,getResult("Inertia", i));
		correlation=Array.concat(correlation,getResult("Correlation", i));
	}
	
	angSecMom_av=newArray(rmc*c);
	invDiffMom_av=newArray(rmc*c);
	contrast_av=newArray(rmc*c);
	energy_av=newArray(rmc*c);
	entropy_av=newArray(rmc*c);
	homogeneity_av=newArray(rmc*c);
	variance_av=newArray(rmc*c);
	shade_av=newArray(rmc*c);
	prominence_av=newArray(rmc*c);
	inertia_av=newArray(rmc*c);
	correlation_av=newArray(rmc*c);

	for(c=1;c<=channels;c++){
		for(i=0; i<rmc;i++){
			angSecMom_sum=0;
			invDiffMom_sum=0;
			contrast_sum=0;
			energy_sum=0;
			entropy_sum=0;
			homogeneity_sum=0;
			variance_sum=0;
			shade_sum=0;
			prominence_sum=0;
			inertia_sum=0;
			correlation_sum=0;
			for(a=0; a<angles.length;a++){
				index=((c-1)*angles.length*rmc)+(a*rmc)+i;
				
				angSecMom_sum=angSecMom_sum+angSecMom[index];
				invDiffMom_sum=invDiffMom_sum+invDiffMom[index];
				contrast_sum=contrast_sum+contrast[index];
				energy_sum=energy_sum+energy[index];
				entropy_sum=entropy_sum+entropy[index];
				homogeneity_sum=homogeneity_sum+homogeneity[index];
				variance_sum=variance_sum+variance[index];
				shade_sum=shade_sum+shade[index];
				prominence_sum=prominence_sum+prominence[index];
				inertia_sum=inertia_sum+inertia[index];
				correlation_sum=correlation_sum+correlation[index];
			}
			index=((c-1)*rmc)+i;
			
			angSecMom_av[index]=angSecMom_sum/angles.length;
			invDiffMom_av[index]=invDiffMom_sum/angles.length;
			contrast_av[index]=contrast_sum/angles.length;
			energy_av[index]=energy_sum/angles.length;
			entropy_av[index]=entropy_sum/angles.length;
			homogeneity_av[index]=homogeneity_sum/angles.length;
			variance_av[index]=variance_sum/angles.length;
			shade_av[index]=shade_sum/angles.length;
			prominence_av[index]=prominence_sum/angles.length;
			inertia_av[index]=inertia_sum/angles.length;
			correlation_av[index]=correlation_sum/angles.length;
		}
	}

	selectWindow("Results"); run("Close");
	IJ.renameResults("Temp","Results");

	for(c=1;c<=channels;c++){
		for(i=0; i<rmc;i++){
			index=(c-1)*rmc+i;
			setResult("AngularSecondMoment", index, angSecMom_av[index]);
			setResult("InverseDifferenceMoment", index, invDiffMom_av[index]);
			setResult("Contrast", index, contrast_av[index]);
			setResult("Energy", index, energy_av[index]);
			setResult("Entropy", index, entropy_av[index]);
			setResult("Homogeneity", index, homogeneity_av[index]);
			setResult("Variance", index, variance_av[index]);
			setResult("Shade", index, shade_av[index]);
			setResult("Prominence", index, prominence_av[index]);
			setResult("Inertia", index, inertia_av[index]);
			setResult("Correlation", index, correlation_av[index]);	
		}
	}
	updateResults();
}


function sortResults()
{
	resultLabels = getResultLabels();
	matrix = results2matrix(resultLabels);
	matrix2results(matrix,resultLabels,channels);
}

function getResultLabels()
{
	selectWindow("Results");
	ls 				= split(getInfo(),'\n');
	rr 				= split(ls[0],'\t'); 
	nparams 		= rr.length-1;			
	resultLabels 	= newArray(nparams);
	for(j=1;j<=nparams;j++){resultLabels[j-1]=rr[j];}
	return resultLabels;
}

function results2matrix(resultLabels)
{
	h = nResults;
	w = resultLabels.length;
	newImage("Matrix", "32-bit Black",w, h, 1);
	matrix = getImageID;
	for(j=0;j<w;j++)
	{
		for(r=0;r<h;r++)
		{
			v = getResult(resultLabels[j],r);
			selectImage(matrix);
			setPixel(j,r,v);
		}
	}
	run("Clear Results");
	return matrix;
}

function matrix2results(matrix,resultLabels,channels)
{
	selectImage(matrix);
	w = getWidth;
	h = getHeight;
	for(c=0;c<channels;c++)
	{
		start = c*h/channels;
		end = c*h/channels+h/channels;
		for(k=0;k<w;k++)
		{
			for(j=start;j<end;j++)
			{
				selectImage(matrix);
				p = getPixel(k,j);
				setResult(resultLabels[k]+"_MC"+c+1,j-start,p); // MC for measurement channel
			}
		}
	}
	selectImage(matrix); close;
	updateResults;
}

function toggleOverlay()
{	
	run("Select None"); 
	roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}

function summarizeResults()
{
	// 	open nuclei results
	run("Results... ", "open=["+nuclei_results+"]");
	nnr 			= nResults;
	nindices		= newArray(nnr);
	resultLabels 	= getResultLabels();
	matrix 			= results2matrix(resultLabels);
	selectWindow("Results"); 
	run("Close");
	for(r=0;r<nnr;r++)
	{
		for(s=0;s<resultLabels.length;s++)
		{
			selectImage(matrix);
			p = getPixel(s,r);
			if(resultLabels[s]!="Cell" && resultLabels[s]!="X_MC1" && resultLabels[s]!="Y_MC1")setResult("Nucl_SC"+nuclei_channel+"_"+resultLabels[s],r,p); // Label all nuclear measured parameters except for the cell or X and Y indices with a "Nucl" prefix
			else if(resultLabels[s]=="X_MC1")setResult("X",r,p);  //exception for X,Y coordinates for ease of tracing-back
			else if(resultLabels[s]=="Y_MC1")setResult("Y",r,p); 
			else setResult(resultLabels[s],r,p);
		}
	}
	updateResults;
	selectImage(matrix); close;
	
	//	append cellular results
	if(File.exists(cells_results))
	{	
		// once in a while a cell index is different from a nuclear index
		for(r=0;r<nnr;r++){nindices[r]=getResult("Cell",r)-1;} 
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+cells_results+"]");
		cnr				= nResults;
		cindices		= newArray(cnr);
		for(r=0;r<cnr;r++){cindices[r]=getResult("Nuclei",r)-1;}
		 
		//Append results based on roi index in cell results 
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(r=0;r<cnr;r++) 
		{
			for(s=0;s<resultLabels.length;s++)
			{
				if(resultLabels[s]!="Nuclei" && resultLabels[s]!="X_MC1" && resultLabels[s]!="Y_MC1")
				{
					selectImage(matrix);
					p = getPixel(s,r);
					setResult("Cell_SC"+cells_channel+"_"+resultLabels[s],cindices[r],p); // Label all cytoplasmic measured parameters with a "Cell" prefix
				}
			}
		}
		updateResults;
		selectImage(matrix); close;
	}
	// 	get distance results (summary currently only implemented for cells or nuclei as secondary objects)
	if (measure_distances)
	{
		if(File.exists(distance_results))
		{
			IJ.renameResults("Results","Temp");
			run("Results... ", "open=["+distance_results+"]");
			min_dist_array = newArray(nResults);
			max_dist_array = newArray(nResults);
			for(r=0;r<nResults;r++)
			{	
				min_dist_array[r] = getResult("Min",r);
				max_dist_array[r] = getResult("Max",r);
			}
			selectWindow("Results"); run("Close");
			IJ.renameResults("Temp","Results");
		}
		else  // no distances because objects were lacking
		{
			min_dist_array = Array.fill(min_dist_array, NaN);
		}
		for(r=0;r<nResults;r++)
		{
			if(secondary_objects =="Nuclei")
			{
				setResult("Nucl_SC"+nuclei_channel+"_MinDist2"+primary_objects,r,min_dist_array[r]);
				setResult("Nucl_SC"+nuclei_channel+"_MaxDist2"+primary_objects,r,max_dist_array[r]);
			}
			if(secondary_objects =="Cells")
			{
				setResult("Cell_SC"+cells_channel+"_MinDist2"+primary_objects,r,min_dist_array[indices[r]]);
				setResult("Cell_SC"+cells_channel+"_MaxDist2"+primary_objects,r,max_dist_array[indices[r]]);
			}
		}
		updateResults;
	}	
	//	append summarized spot results
	if(File.exists(spots_a_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+spots_a_results+"]");
		// if distance measurement was done for spots add them here
		if (measure_distances && secondary_objects=="Spots_a")
		{
			for(r=0;r<min_dist_array.length;r++)
			{
				setResult("MinDist2"+primary_objects,r,min_dist_array[r]);
				setResult("MaxDist2"+primary_objects,r,max_dist_array[r]);
			}
		}
		updateResults;
		snr 			= nResults;
		nindices 		= newArray(snr);
		cindices 		= newArray(snr);
		for(j=0;j<snr;j++)
		{
			nindices[j] = getResult("Nucleus",j)-1;
			cindices[j] = getResult("Cell",j)-1;
		}	
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(s=0;s<resultLabels.length;s++)
		{
			if(resultLabels[s] != "Nucleus" && resultLabels[s] != "Cell")
			{
				nvalues 	= newArray(nnr);
				cvalues 	= newArray(nnr);
				nnumber 	= newArray(nnr);
				cnumber 	= newArray(nnr);
				for(r=0;r<snr;r++)
				{
					selectImage(matrix);
					p = getPixel(s,r);
					if(nindices[r]>=0)
					{
						nvalues[nindices[r]] += p;  
						nnumber[nindices[r]] += 1;	
					}
					if(cindices[r]>=0)
					{
						cvalues[cindices[r]] += p;  
						cnumber[cindices[r]] += 1;	
					}
				}
				
				for(r=0;r<nnr;r++)
				{
					setResult("Spot_SC"+spots_a_channel+"_NrPerNuc",r,nnumber[r]);
					setResult("Spot_SC"+spots_a_channel+"_"+resultLabels[s]+"_SumPerNuc",r,nvalues[r]);              
					setResult("Spot_SC"+spots_a_channel+"_"+resultLabels[s]+"_MeanPerNuc",r,nvalues[r]/nnumber[r]);
					if(segment_cells)
					{
						setResult("Spot_SC"+spots_a_channel+"_NrPerCell",r,cnumber[r]);
						setResult("Spot_SC"+spots_a_channel+"_"+resultLabels[s]+"_SumPerCell",r,cvalues[r]);
						setResult("Spot_SC"+spots_a_channel+"_"+resultLabels[s]+"_MeanPerCell",r,cvalues[r]/cnumber[r]);
					}
				}
			}
		}
		selectImage(matrix); close;
		updateResults();
	}
	if(File.exists(spots_b_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+spots_b_results+"]");
		// if distance measurement was done for spots add them here
		if (measure_distances && secondary_objects=="Spots_b")
		{
			for(r=0;r<min_dist_array.length;r++)
			{
				setResult("MinDist2"+primary_objects,r,min_dist_array[r]);
				setResult("MaxDist2"+primary_objects,r,max_dist_array[r]);
			}
		}
		updateResults;
		snr 			= nResults;
		nindices 		= newArray(snr);
		cindices 		= newArray(snr);
		for(j=0;j<snr;j++)
		{
			nindices[j] = getResult("Nucleus",j)-1;
			cindices[j] = getResult("Cell",j)-1;
		}	
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(s=0;s<resultLabels.length;s++)
		{
			if(resultLabels[s] != "Nucleus" && resultLabels[s] != "Cell")
			{
				nvalues 	= newArray(nnr);
				cvalues 	= newArray(nnr);
				nnumber 	= newArray(nnr);
				cnumber 	= newArray(nnr);
				for(r=0;r<snr;r++)
				{
					selectImage(matrix);
					p = getPixel(s,r);
					if(nindices[r]>=0)
					{
						nvalues[nindices[r]] += p;  
						nnumber[nindices[r]] += 1;	
					}
					if(cindices[r]>=0)
					{
						cvalues[cindices[r]] += p;  
						cnumber[cindices[r]] += 1;	
					}
				}
				for(r=0;r<nnr;r++)
				{
					setResult("Spot_SC"+spots_b_channel+"_NrPerNuc",r,nnumber[r]);
					setResult("Spot_SC"+spots_b_channel+"_"+resultLabels[s]+"_SumPerNuc",r,nvalues[r]);              
					setResult("Spot_SC"+spots_b_channel+"_"+resultLabels[s]+"_MeanPerNuc",r,nvalues[r]/nnumber[r]);
					if(segment_cells)
					{
						setResult("Spot_SC"+spots_b_channel+"_NrPerCell",r,cnumber[indices[r]]);
						setResult("Spot_SC"+spots_b_channel+"_"+resultLabels[s]+"_SumPerCell",r,cvalues[r]);
						setResult("Spot_SC"+spots_b_channel+"_"+resultLabels[s]+"_MeanPerCell",r,cvalues[r]/cnumber[r]);
					}
				}
			}
		}
		selectImage(matrix); close;
		updateResults();
	}
	// fold results
	if(File.exists(folds_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+folds_results+"]");
		fnr 			= nResults;
		nindices 		= newArray(fnr);
		cindices 		= newArray(fnr);
		for(j=0;j<fnr;j++)
		{
			nindices[j] = getResult("Nucleus",j)-1;
			cindices[j] = getResult("Cell",j)-1;
		}	
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(s=0;s<resultLabels.length;s++)
		{
			if(resultLabels[s] != "Nucleus" && resultLabels[s] != "Cell")
			{
				nvalues 	= newArray(nnr);
				cvalues 	= newArray(nnr);
				nnumber 	= newArray(nnr);
				cnumber 	= newArray(nnr);
				for(r=0;r<fnr;r++)
				{
					selectImage(matrix);
					p = getPixel(s,r);
					if(nindices[r]>=0)
					{
						nvalues[nindices[r]] += p;  
						nnumber[nindices[r]] += 1;	
					}
					if(cindices[r]>=0)
					{
						cvalues[cindices[r]] += p;  
						cnumber[cindices[r]] += 1;	
					}
				}
				for(r=0;r<nnr;r++)
				{
					setResult("Fold_SC"+folds_channel+"_NrPerNuc",r,nnumber[r]);
					setResult("Fold_SC"+folds_channel+"_"+resultLabels[s]+"_SumPerNuc",r,nvalues[r]);              
					setResult("Fold_SC"+folds_channel+"_"+resultLabels[s]+"_MeanPerNuc",r,nvalues[r]/nnumber[r]);
					if(segment_cells)
					{
						setResult("Fold_SC"+folds_channel+"_NrPerCell",r,cnumber[r]);
						setResult("Fold_SC"+folds_channel+"_"+resultLabels[s]+"_SumPerCell",r,cvalues[r]);
						setResult("Fold_SC"+folds_channel+"_"+resultLabels[s]+"_MeanPerCell",r,cvalues[r]/cnumber[r]);
					}
				}
			}
		}
		selectImage(matrix); close;
		updateResults();
	}
	selectWindow("Results"); saveAs("Measurements",results);
}
	
function createOverlay(names)
{
	setForegroundColor(25, 25, 25);
	fields = names.length;
	print(fields,"images");
	for(i=0;i<fields;i++)
	{
		prefix = names[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id = getImageID;
		Stack.getDimensions(w,h,channels,slices,frames); 
		if(!Stack.isHyperStack && channels == 1)
		{
			channels = slices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
		}
		id = getImageID;

		if(segment_nuclei)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(nuclei_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",nuclei_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
	
		if(segment_cells)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(cells_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",cells_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}

		if(segment_spots && spots_a_channel>0)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_a_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_a_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
		
		if(segment_spots && spots_b_channel>0)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_b_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_b_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
		
		if(segment_folds)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(folds_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",folds_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}

		if(segment_plaques)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(plaques_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",plaques_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
	}
	run("Concatenate...", "all_open title=[Concatenated Stacks]");
	Stack.getDimensions(w,h,newchannels,slices,frames);
	for(c=1;c<=channels;c++){Stack.setChannel(c);Stack.setFrame(round(frames/2));resetMinAndMax;}
	range = pow(2,bitDepth);
	for(c=channels+1;c<=newchannels;c++){Stack.setChannel(c);setMinAndMax(0,range/2);}
	run("Make Composite");
}