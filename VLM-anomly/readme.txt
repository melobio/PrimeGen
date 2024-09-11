# Environment requirement: vllm, flash-attention, swift, CVAT, deepspeed, pandas, openai

# Special requirement: OpenAI key is needed for data augmentation

# Experts use CVAT for data labeling
# After labeling, use `extract_annotation.py` for processing  
# After processing, use `gen_qa.py` to generate VQA (Visual Question Answering) data  
# Once VQA data is generated, use `make_swift_data.py` to divide the data  
# After data division, use `swift_phi3_sft_lora.sh` for LoRA training  
# After training, use `swift_phi3_merging.sh` to merge the model weights  
# After merging, use `real_infer_phi3_vllm.py` for inference
