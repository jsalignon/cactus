

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [Data](/docs/2_Install/Data.md), [Test_datasets](/docs/2_Install/Test_datasets.md)
* [Run](/docs/3_Run/3_Run.md): [Input Data](/docs/3_Run/Input_data.md), [Input Files](/docs/3_Run/Input_files.md), [Parameters](/docs/3_Run/Parameters.md)


[](END_OF_MENU)






Cactus was built with the goal to make the installation as pain-free as possible. This was achieved by using only tools within containers for all analysis. Thereby, installation of most tools is done simply by downloading containers. However, two key dependencies are still necessary. These are Nextflow (the pipeline tool) and SingularityCE (the container tool).

In case of issue with the pipeline and/or for reproducibility purporsess, it might be useful to try to run Cactus with the version of these tools that was used when developing the current Cactus release. These are:
  - SingularityCE: singularity-ce version 3.10.0+91-g13f189977 (released on May 17, 2022)
  - Nextflow: version 22.05.0-edge build 5704 (released on May 25, 2022 )
