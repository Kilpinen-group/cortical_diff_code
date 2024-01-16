
#!/bin/sh

cd cranger/Batch3-4_091221-161221/multi_220121_A00464_0441_BHVYHGDSX2
export MPSTKZ=8M
cellranger-6.1.2/cellranger multi --id=Day_20_Pool1 --csv=Multi_Config_Day_20_Pool1.csv
cellranger-6.1.2/cellranger multi --id=Day_20_Pool2 --csv=Multi_Config_Day_20_Pool2.csv


