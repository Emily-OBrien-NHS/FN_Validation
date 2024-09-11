Code to create pathway definition and other outputs from events file and perform some validation.

main.py - The main file to run with inputted events data to create the pathway definition and other input files using various filtering/settings.
	  Option to include diagnostic, obs and admission data.
data_cleaning_and_transformation.py - The file containing all the data cleaning/transformation functions.  Some initial cleaning functions at the start,
				      with the main_clense_and_transform_data function at the end, combining most of these functions.
pathway_definitions.py - The file containing all the functions to create the pathway definitions and directed flow diagrams, with
			 main_generate_dfg_and_pathway_definitions function at the end, combining these functions.
process_durations.pt - The file containing all the functions to calculate process durations and fit log normal distributions, with
 		       main_generate_histogram_and_process_durations function at the end, combining all these functions.
config.py - The file containing all the input/config settings to run this code.

validation.py - The file to perform validation on an inputted full scenario, to find any processes/resources/locations that may not match up. 
		Creates an outputted text file of all the possible mismatches/missing inputs.
