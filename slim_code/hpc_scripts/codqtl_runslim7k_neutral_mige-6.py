###Python/pyslim code for running cod polygenic adaptation + evolutionary rescue simulations 

##Import needed packages
import subprocess
import multiprocessing
import os
import sys

##set up needed common paths
dir="/scratch/br450/codslim/"

slim_template= dir + "codqtl_reduced_neutral_slo_template.slim"

output_dir= dir + "K7000_m1e-6_neutral/"

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

def run_sim(i):
	Lof07out=output_dir+"i"+str(i)+"_Lof07.vcf"
	Can40out=output_dir+"i"+str(i)+"_Can40.vcf"
	Lof11out=output_dir+"i"+str(i)+"_Lof11.vcf"
	Lof14out=output_dir+"i"+str(i)+"_Lof14.vcf"
	Can13out=output_dir+"i"+str(i)+"_Can13.vcf"
	simout=output_dir+"i"+str(i)+"_simout.txt"
	migout=output_dir+"i"+str(i)+"_migout.txt"
	try:
		subprocess.check_output(["/home/br450/build_082622/slim","-m","-d","K=7000","-d","m=0.000001","-d","simLof07file='"+Lof07out+"'","-d","simCan40file='"+Can40out+"'","-d","simLof11file='"+Lof11out+"'","-d","simLof14file='"+Lof14out+"'","-d","simCan13file='"+Can13out+"'","-d","simoutfile='"+simout+"'","-d","migoutfile='"+migout+"'",slim_template])
	except Exception as e:
		return e


##run twenty iterations
a_pool = multiprocessing.Pool()
result = a_pool.map(run_sim, range(1,21))
