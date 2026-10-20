import sys, types, runpy
fai, bed, out_bed, out_pri, script = sys.argv[1:6]
snakemake = types.SimpleNamespace(input=[fai, bed], output=[out_bed, out_pri])
g = runpy.run_path(script, init_globals={"snakemake": snakemake}, run_name="__main__")
