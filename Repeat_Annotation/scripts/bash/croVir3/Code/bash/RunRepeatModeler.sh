#!/bin/bash

: <<'ScriptDescription'

This script is designed to run RepeatModeler. I am following the method found in Daren Card's blog post:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

ScriptDescription

# Ensure that the script returns the exit code of a failed command in a pipeline
set -o pipefail

# Activate the mamba environment containing RepeatModeler
eval "$(micromamba shell hook --shell bash)"
micromamba activate repeatmask

# Set the start time of the script
start_time=$(date '+%Y-%m-%d %H:%M:%S')

# Set ntfy.sh topic for notifications
ntfy_topic="kaas-ballard-Klauber-scripts-27857274017852061578"

# NOTE: Change this if needed
# Set the number of threads for RepeatModeler
threads=10

# SET: Change this each time you run this script
# Set the name for the RepeatModeler database
species_name="Crotalus_viridis_viridis"

# SET: Change this each time you run this script
# Set the family name for the species
family="Viperidae"

# SET: Change this each time you run this script
# Set the genus name for the species
genus="Crotalus"

# SET: Change this each time you run the script
# Reference genome
reference_genome="$HOME/Documents/Kaas/SquamateAlignments/Data/Reference_Genomes/CastoeGenomes/Crotalus/Crotalus_viridis_viridis/Assembly/Cvv_2017_genome_with_myos.fasta"

# NOTE: Change this each time you run the script
# Set the output directory
output_directory="$HOME/Documents/Kaas/SquamateAlignments/SoftMasking/$family/$genus/$species_name/croVir3/Results/0_RepeatModeler"

# Create log directory under the output directory if it does not exist
[ ! -d "$output_directory/Logs" ] && mkdir -p "$output_directory/Logs"

# Send myself a notification that the script is starting
curl -d "🔔 Starting RepeatModeler for $species_name at $(date). Reference genome: $reference_genome. Check logs at $output_directory/Logs/ for details." \
	ntfy.sh/"$ntfy_topic"

# Change the directory to the output directory so that the RM_ files are created there
cd "$output_directory" || { echo "Failed to change directory to $output_directory"; exit 1; }

# Step #1: Build a new RepeatModeler database for the reference genome
# Note: BuildDatabase can build directories that don't exist yet
if [ ! -s "$output_directory/$species_name.nsq" ]; then
	echo "Warning: Database at the specified path does not exist or is empty. Building a new database now."
	BuildDatabase -name "$output_directory/$species_name" \
		"$reference_genome" 2>&1 | tee "$output_directory/Logs/BuildDatabase.log"
else
	echo "Database already exists at the specified path and is not empty. Skipping."
fi

# Step #2: Detect if the database was created and then run RepeatModeler
# Create a variable to hold the exit code of RepeatModeler
repeat_modeler_exit_code=0
if [ ! -s "$output_directory/$species_name.nsq" ]; then
	echo "Error: Database was not created or is empty. Exiting."
	exit 1
else
	# Find the most recent recovery directory, if any
	recover_dir=$(ls -td "$output_directory"/RM_*/ 2>/dev/null | head -n 1)
	if [ -n "$recover_dir" ]; then
		echo "Previous RepeatModeler run found. Attempting to recover from $recover_dir."

		# Run RepeatModeler
		RepeatModeler \
			-threads "$threads" \
			-database "$output_directory/$species_name" \
			-recoverDir "$recover_dir" 2>&1 | tee "$output_directory/Logs/RepeatModeler.log"
		repeat_modeler_exit_code=$?
	else
		echo "No previous RepeatModeler run found. Starting a new run."

		# Run RepeatModeler
		RepeatModeler \
			-threads "$threads" \
			-database "$output_directory/$species_name" 2>&1 | tee "$output_directory/Logs/RepeatModeler.log"
		repeat_modeler_exit_code=$?
	fi
fi

# Step #3: Cleanup temporary files in the RM_*** directory
if [ $repeat_modeler_exit_code -eq 0 ]; then
	echo "RepeatModeler completed successfully."
	# Add a safety check to ensure the output directory is not something dangerous
	if [[ -n "$output_directory" && "$output_directory" != "/" && "$output_directory" != "$HOME" && "$output_directory" != "." && "$output_directory" != ".." ]]; then
		# RepeatModeler completed successfully, cleanup any extraneous RM_ directories
		echo "Cleaning up temporary RM_* directories..."
		for dir in "$output_directory"/RM_*; do
			# This check handles the case where no directories match the glob
			if [ -d "$dir" ]; then
				echo "Removing $dir"
				rm -rf "$dir"
			fi
		done
	fi
else
	echo "RepeatModeler did not complete successfully, exiting and sending notification to user."
	curl -d "❌ FAILED: RepeatModeler failed for $species_name at $(date) prior to cleanup of temporary files. Check logs for errors." \
		ntfy.sh/"$ntfy_topic"
	exit 1
fi

# Set the end time of the script
end_time=$(date '+%Y-%m-%d %H:%M:%S')

# Calculate the duration of the script
runtime=$(date -u -d "$end_time" +"%s")-$(( $(date -u -d "$start_time" +"%s") ))

# Tell the user in the log file and term how long the script took to run
echo "Total runtime: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) seconds"

# Send a notification that the script has finished
if [ $repeat_modeler_exit_code -eq 0 ] && \
	[ -s "$output_directory/$species_name.nsq" ] && \
	[ -s "$output_directory/$species_name-families.fa" ] && \
	[ -s "$output_directory/$species_name-families.stk" ] && \
	[ -s "$output_directory/$species_name.njs" ]; then # if [ $? -eq 0 ] checks if the last command was successful
	curl -d "✅ SUCCESS: RepeatModeler completed for $species_name at $(date). Total runtime: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) secondsCheck logs at $output_directory/Logs/" \
		ntfy.sh/"$ntfy_topic"
else
	curl -d "❌ FAILED: RepeatModeler failed for $species_name at $(date). Check logs for errors." \
		ntfy.sh/"$ntfy_topic"
fi