#The user selects the number of nodes and device ID based on the configuration, while also providing the Base address and the address for data construction
NPROC_PER_NODE=4 CUDA_VISIBLE_DEVICES=1,2,3,4 swift sft  \
  --model_type phi3-vision-128k-instruct  \
  --model_id_or_path "base_path" \
  --custom_train_dataset_path train.jsonl \
  --custom_val_dataset_path val.jsonl \
  --batch_size 1 \
  --num_train_epochs 10 \
  --ddp_find_unused_parameters true \
  --save_total_limit -1 \
  --save_strategy epoch \
  --use_flash_attn true \
  --safe_serialization false \
  --evaluation_strategy epoch
