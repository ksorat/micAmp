#module swap intel/15.0.0
module swap intel/14.0.2
#module unload mvapich2/2.0/INTEL-14.0.2
#module load impi

export OMP_NUM_THREADS=16
export KMP_AFFINITY=scatter
export MIC_ENV_PREFIX=MIC
export MIC_OMP_NUM_THREADS=240
export MIC_KMP_AFFINITY=granularity=fine,scatter
export OFFLOAD_INIT=on_start
export OFFLOAD_REPORT=0
export KMP_STACKSIZE=128M
export MIC_KMP_STACKSIZE=32M
export MIC_STACKSIZE=256M
export MIC_USE_2MB_BUFFERS=512k