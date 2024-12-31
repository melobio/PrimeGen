## Setting up the Environment

To set up the development environment for this project, we use `conda`. The environment configuration is saved in an `environment.yml` file, which includes both `conda` and `pip` dependencies.

### Prerequisites

Make sure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) installed on your machine.

### 1. Create the Conda Environment

Run the following command to create the environment using the `environment.yml` file:

```bash
conda env create -f environment.yml


# Special requirement: OpenAI key is needed for data augmentation

1) Experts use CVAT for data labeling.
2) Use `data_generation.py` for data generation.
3) Use `create_stage1_data.py` for generating stage 1 data.
4) Use `create_stage2_data.py` for generating stage 2 data.
5) Use `qwen_train_stage1.sh` for training the stage 1 model.
6) Use `merging_lora.sh` for merging the stage 1 model with LoRA.
7) For the base stage LoRA model, use `qwen_train_stage2.sh` for stage 2 training.
8) Use `merging_lora.sh` again for merging the stage 2 model with LoRA.
9) Use `deploy.sh` to deploy the stage 1 model.
10) Use `deploy_stage2.sh` to deploy the stage 2 model.
11) Use `infer_stage1.py` for inference on the stage 1 test dataset.
12) Use `infer_stage2.py` for inference on the stage 2 test dataset.