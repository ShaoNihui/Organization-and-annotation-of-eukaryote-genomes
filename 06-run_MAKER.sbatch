#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=120G
#SBATCH -p pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=MAKER_N13

#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/MAKER_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/MAKER_%j.err

set -euo pipefail

# 课程目录（容器和 RepeatMasker 在这里）
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

# 你的工作目录（MAKER 控制文件和输出）
WORKDIR="/data/users/nshao/organization_annotation_course/gene_annotation"

# 你的整个用户数据根目录 —— 里面包含：
# - organization_annotation_course
# - assembly_annotation_course
# - 以及其他个人项目
USERROOT="/data/users/nshao"

REPEATMASKER_DIR="$COURSEDIR/softwares/RepeatMasker"
export PATH="$PATH:$REPEATMASKER_DIR"

cd "$WORKDIR"

echo "[INFO] Starting MAKER MPI run in $WORKDIR on node $(hostname) at $(date)"

module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

mpiexec --oversubscribe -n 50 apptainer exec \
  --bind $SCRATCH:/TMP \
  --bind "$USERROOT" \
  --bind "$COURSEDIR" \
  --bind "$AUGUSTUS_CONFIG_PATH" \
  --bind "$REPEATMASKER_DIR" \
  "$COURSEDIR/containers/MAKER_3.01.03.sif" \
  maker -mpi -ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "[INFO] MAKER finished at $(date)"