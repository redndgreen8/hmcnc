import os
import tempfile
import subprocess
import os.path
import json
import argparse


#parser = argparse.ArgumentParser(description="Location of boost include folder")
#parser.add_argument("--boost",help="Location of boost install/include folder. \
#	Most likely {anaconda install}/envs/\{proj_env\}/include")

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)


boost=config['boost']


rule all:
	input:
		samtobed="samToBed",
		vit="viterbi",

rule stb:
	input:
		sTb="SamToBed.cpp",
	output:
		samtobed="samToBed",
	shell:"""

		g++ -02 {input.sTb} -o {output.samtobed}
"""


rule vit:
	input:
		vitt="viterbi.noclip.cpp",
	output:
		samtobed="viterbi",
	params:
		proj_env=boost,

	shell:"""
g++ -W -I {proj_env}/include {input.vitt} -o {samtobed}
"""


	