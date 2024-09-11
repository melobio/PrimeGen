#The user provides the addresses of the training models whose weights need to be merged.
CUDA_VISIBLE_DEVICES=0 swift export \
	--ckpt_dir "lora_path" \
	--merge_lora true \
	--safe_serialization false
