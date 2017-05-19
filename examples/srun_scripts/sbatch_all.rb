srun_dir = "/home/agossman/github/FDRcorrectedSCCA/examples/srun_scripts/"
srun_files = Dir[srun_dir + "/*.srun"]

srun_files.each do |file_name|
  system("sbatch #{file_name}")
end
