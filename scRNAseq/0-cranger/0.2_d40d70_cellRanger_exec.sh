#!/bin/sh

cd cranger/Batch1-2_220921-221021/count_211123_A00464_0418_BHNVTFDSX2
qsub -t 1-5 d40d70_cellRanger.sh d40d70_cellRanger.txt
